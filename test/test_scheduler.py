"""Slurm contract tests; subprocesses are mocked, never real submissions."""
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from py4siesta import scheduler, tool_cli
from py4siesta.operations import JobSubmissionOperation
from py4siesta_agent.scheduler import SchedulerManager


def reply(stdout='', stderr=''):
    return SimpleNamespace(stdout=stdout, stderr=stderr)


class SchedulerTests(unittest.TestCase):
    def test_explicit_script_is_relative_to_case_and_needs_no_origin_or_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp) / 'case with spaces'
            case.mkdir()
            script = case / 'run script.sh'
            script.write_text('#!/bin/sh\ntrue\n')
            with mock.patch.object(scheduler.subprocess, 'run', return_value=reply('123;cluster-a\n')) as run:
                result = scheduler.submit_job(str(case), script.name)
            self.assertEqual(result['job_id'], '123')
            self.assertEqual(result['cluster'], 'cluster-a')
            self.assertEqual(result['status'], 'submitted')
            self.assertEqual(run.call_args.args[0], ['sbatch', '--parsable', str(script)])
            self.assertEqual(run.call_args.kwargs['cwd'], str(case))
            self.assertFalse((case / 'origin').exists())

    def test_missing_script_does_not_submit(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(scheduler.subprocess, 'run') as run:
            with self.assertRaises(FileNotFoundError):
                scheduler.submit_job(tmp, 'missing.sh')
            run.assert_not_called()

    def test_ambiguous_submission_preserves_evidence_in_cli_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'run.sh').write_text('#!/bin/sh\ntrue\n')
            with mock.patch.object(scheduler.subprocess, 'run', return_value=reply('unexpected response\n')) as run:
                result = tool_cli.execute(['job-submit', '--case', tmp, '--script', 'run.sh'])
            self.assertFalse(result['ok'])
            details = result['error']['details']
            self.assertEqual(details['status'], 'submission_unknown')
            self.assertFalse(details['retry_safe'])
            self.assertEqual(details['stdout'], 'unexpected response\n')
            self.assertEqual(run.call_count, 1)

    def test_submission_timeout_is_uncertain_and_never_retried(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'run.sh').write_text('#!/bin/sh\ntrue\n')
            error = subprocess.TimeoutExpired(['sbatch'], 30, output=b'partial', stderr=b'connection lost')
            with mock.patch.object(scheduler.subprocess, 'run', side_effect=error) as run:
                with self.assertRaises(scheduler.SchedulerError) as caught:
                    scheduler.submit_job(tmp, 'run.sh')
            self.assertEqual(caught.exception.details['stdout'], 'partial')
            self.assertEqual(caught.exception.details['status'], 'submission_unknown')
            self.assertEqual(run.call_count, 1)

    def test_active_status_is_read_only_and_preserves_raw_state(self):
        with mock.patch.object(scheduler.subprocess, 'run', return_value=reply('CLUSTER: cluster-a\n123|RUNNING\n')) as run:
            result = scheduler.job_status('123', 'cluster-a')
        self.assertEqual(result['status'], 'running')
        self.assertEqual(result['raw_state'], 'RUNNING')
        self.assertIsNone(result['exit_code'])
        self.assertEqual(run.call_count, 1)
        self.assertEqual(run.call_args.args[0], ['squeue', '-h', '-j', '123', '-o', '%i|%T', '--clusters', 'cluster-a'])

    def test_accounting_selects_allocation_not_batch_step(self):
        with mock.patch.object(scheduler.subprocess, 'run', side_effect=[
            reply(), reply('123.batch|FAILED|1:0\n123|COMPLETED|0:0\n')
        ]) as run:
            result = scheduler.job_status('123')
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['exit_code'], '0:0')
        self.assertEqual(result['source'], 'sacct')
        self.assertEqual([c.args[0][0] for c in run.call_args_list], ['squeue', 'sacct'])

    def test_live_queue_rejection_can_still_resolve_completed_job(self):
        with mock.patch.object(scheduler.subprocess, 'run', side_effect=[
            subprocess.CalledProcessError(1, ['squeue'], stderr='Invalid job id specified'),
            reply('123|COMPLETED|0:0\n')
        ]):
            self.assertEqual(scheduler.job_status('123')['status'], 'completed')

    def test_array_task_and_cancelled_suffix(self):
        with mock.patch.object(scheduler.subprocess, 'run', side_effect=[reply(), reply('123_4|CANCELLED by 1045|0:15\n')]):
            result = scheduler.job_status('123_4')
        self.assertEqual(result['status'], 'cancelled')
        self.assertEqual(result['raw_state'], 'CANCELLED by 1045')
        self.assertEqual(result['exit_code'], '0:15')

    def test_unknown_records_and_states_are_not_failed_or_completed(self):
        for accounting in ['', '123|FUTURE_STATE|\n', '123|FAILED|1:0\n123|COMPLETED|0:0\n']:
            with self.subTest(accounting=accounting), mock.patch.object(scheduler.subprocess, 'run', side_effect=[reply(), reply(accounting)]):
                self.assertEqual(scheduler.job_status('123')['status'], 'unknown')
        with mock.patch.object(scheduler.subprocess, 'run', return_value=reply('123_1|RUNNING\n123_2|PENDING\n')):
            self.assertEqual(scheduler.job_status('123')['status'], 'unknown')

    def test_query_error_does_not_report_job_failure(self):
        error = subprocess.CalledProcessError(1, ['squeue'], stderr='controller unavailable')
        with mock.patch.object(scheduler.subprocess, 'run', side_effect=error) as run:
            result = tool_cli.execute(['job-status', '--job-id', '123'])
        self.assertFalse(result['ok'])
        self.assertEqual(result['error']['details']['status'], 'query_error')
        self.assertEqual(result['error']['details']['stderr'], 'controller unavailable')
        self.assertEqual(run.call_count, 2)

    def test_cancel_requests_exact_job_without_claiming_completion(self):
        with mock.patch.object(scheduler.subprocess, 'run', return_value=reply()) as run:
            result = scheduler.cancel_job('123_4', 'cluster-a')
        self.assertEqual(result['status'], 'cancel_requested')
        self.assertEqual(run.call_args.args[0], ['scancel', '--clusters', 'cluster-a', '123_4'])

    def test_invalid_selectors_never_reach_slurm(self):
        with mock.patch.object(scheduler.subprocess, 'run') as run:
            for bad in ['--all', '1,2', '0', '-1', '1;rm', '1.batch']:
                with self.subTest(bad=bad), self.assertRaises(ValueError):
                    scheduler.cancel_job(bad)
            with self.assertRaises(ValueError):
                scheduler.job_status('123', '--all')
            run.assert_not_called()

    def test_legacy_submission_retains_first_directory_and_all_script_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for base in ['01.a', '01.b']:
                for case in ['case', 'optimized_structure']:
                    directory = root / base / case
                    directory.mkdir(parents=True)
                    for script in ['slm_a', 'slm_b', 'ignored.sh']:
                        (directory / script).write_text('#!/bin/sh\n')
            calls = []
            def fake(command, **kwargs):
                calls.append((Path.cwd(), command, kwargs))
            with mock.patch.object(scheduler.subprocess, 'run', side_effect=fake):
                JobSubmissionOperation(SimpleNamespace(root=root)).run('kpt')
            self.assertEqual(len(calls), 4)
            self.assertTrue(all(directory.parent == root/'01.a' for directory, _, _ in calls))
            self.assertTrue(all(command[0] == 'sbatch' and '--parsable' not in command and kwargs == {'check': True}
                                for _, command, kwargs in calls))

    def test_agent_imports_shared_backend_and_retains_budget_on_query_error(self):
        from py4siesta_agent.scheduler import SlurmBackend, SchedulerError
        self.assertIs(SlurmBackend, scheduler.SlurmBackend)
        self.assertIs(SchedulerError, scheduler.SchedulerError)
        backend = mock.Mock(name='backend')
        backend.status.side_effect = scheduler.SchedulerError('controller unavailable')
        manager = SchedulerManager(backend=backend, max_total_nodes=1)
        active = dict(status='running', requested_nodes=1, scheduler_job_id='123')
        pending = dict(status='pending', requested_nodes=1)
        manager.update([active, pending])
        self.assertEqual(active['status'], 'running')
        self.assertIn('query_error', active)
        self.assertNotIn('completed_at', active)
        backend.submit.assert_not_called()
        backend.status.side_effect = None
        backend.status.return_value = 'running'
        manager.update([active])
        self.assertNotIn('query_error', active)


if __name__ == '__main__':
    unittest.main()
