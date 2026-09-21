import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from nanocore import siesta
from py4siesta import post_process


class BandTests(unittest.TestCase):
    def test_legacy_spin_blocks_and_multiline_input(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'sample.bands'
            for nspin in (1, 2):
                nbands = 12
                values = np.arange(nbands * nspin, dtype=float)
                lines = ['5.5', '0 1', '0 24', '%d %d 2' % (nbands, nspin)]
                for k in (0, 1):
                    for start in range(0, len(values), 10):
                        prefix = '%d ' % k if start == 0 else ''
                        lines.append(prefix + ' '.join(str(v + k) for v in values[start:start + 10]))
                lines.extend(['2', "0 'G'", "1 'X'"])
                path.write_text('\n'.join(lines) + '\n')
                with mock.patch.object(siesta.os, 'system', side_effect=AssertionError('external utility')):
                    result = siesta.get_band(None, None, label=str(path.with_suffix('')))
                    data = siesta.get_band(None, None, bands_path=path, return_data=True)
                self.assertEqual(len(result), 2 * nspin)
                for spin in range(nspin):
                    np.testing.assert_allclose(result[2 * spin], [[0, 1]] * nbands)
                    expected = np.array([values + k - 5.5 for k in (0, 1)]).T
                    np.testing.assert_allclose(result[2 * spin + 1], expected[spin * nbands:(spin + 1) * nbands])
                self.assertEqual(data['vbm'], 5.0)
                self.assertEqual(data['bandgap'], 1.0)
                self.assertEqual(data['labels'], ['G', 'X'])
                self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_wrapper_delegates_and_formats_labels(self):
        data = dict(kpath=np.array([0]), energies=np.array([[[1]]]), nbands=1, nkpoints=1,
                    nspin=1, special_k=np.array([0]), labels=['G'],
                    fermi_level=0, bandgap=0, vbm=0)
        with tempfile.NamedTemporaryFile(suffix='.bands') as source:
            with mock.patch.object(siesta, 'get_band', return_value=data) as reader:
                result = post_process.read_band_structure(source.name)
            reader.assert_called_once_with(file_path=Path(source.name), return_data=True)
        self.assertIsInstance(result, post_process.BandStructureData)
        self.assertEqual(result.labels, [r'$\Gamma$'])

    def test_rerun_preserves_input_and_post_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                Path('path.fdf').write_text('BandLinesScale pi/a\n')
                Path('RUN.fdf').write_text('Existing input\n')
                simulation = mock.Mock()
                simulation._params = {'Label': 'siesta'}
                with mock.patch.object(siesta, '_read_band_structure', return_value={}) as reader:
                    self.assertEqual(siesta.get_band(simulation, 'path.fdf', rerun=1, return_data=True), {})
                simulation.run.assert_called_once_with(mode='POST')
                reader.assert_called_once_with(Path('siesta.bands'))
                self.assertEqual(Path('RUN.fdf').read_text(), 'Existing input\nBandLinesScale pi/a\nWriteBands            T')
            finally:
                os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
