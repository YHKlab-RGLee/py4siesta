"""User-invoked registration; never called as part of calculation workflows."""
import argparse
import copy
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true', help='show changes without writing')
    args = parser.parse_args()
    try:
        import tomllib as toml
    except ImportError:
        try:
            import tomli as toml
        except ImportError:
            try:
                import toml
            except ImportError:
                raise ValueError('Use Python 3.11+, or install tomli/toml in the selected Python environment.')

    source = Path(__file__).resolve().parents[1] / 'py4siesta'
    checkout = source.parents[2]
    target = Path.home() / '.agents/skills/py4siesta'
    state = Path(os.environ.get('PY4SIESTA_SKILL_STATE') or
                 str(Path.home() / '.local/state/py4siesta-skill')).expanduser().resolve()
    config = Path(os.environ.get('CODEX_HOME') or str(Path.home() / '.codex')).expanduser() / 'config.toml'
    if not (source / 'SKILL.md').is_file():
        raise ValueError('Bundled SKILL.md is missing.')
    if state == checkout or checkout in state.parents or state == target.resolve() or target.resolve() in state.parents:
        raise ValueError('Personal state must be outside the repository and installed skill tree.')
    if os.path.lexists(str(target)) and not (target.is_symlink() and target.resolve() == source):
        raise ValueError('Registration destination already exists: {}. Move it manually before registering.'.format(target))
    if config.is_symlink():
        raise ValueError('Config is a symlink; update its target manually: {}'.format(config))
    original = config.read_text(encoding='utf-8') if config.exists() else ''
    data = toml.loads(original)
    expected = copy.deepcopy(data)
    section = expected.setdefault('sandbox_workspace_write', {})
    roots = section.setdefault('writable_roots', [])
    if not isinstance(roots, list) or not all(isinstance(p, str) for p in roots):
        raise ValueError('Existing writable_roots must be an array of strings.')
    updated = original
    if str(state) not in roots:
        roots.append(str(state))
        value = json.dumps(roots, ensure_ascii=True)
        # Accept an edit only if parsing proves that all other settings survive.
        candidates = []
        headers = list(re.finditer(r'^\s*\[sandbox_workspace_write\][ \t]*(?:#[^\n]*)?$', original, re.M))
        if headers:
            header = headers[0]
            tail = original[header.end():]
            next_header = re.search(r'^\s*\[', tail, re.M)
            end = header.end() + next_header.start() if next_header else len(original)
            block = original[header.end():end]
            key = re.search(r'^[ \t]*writable_roots[ \t]*=', block, re.M)
            if key:
                start = header.end() + key.start()
                for line_end in re.finditer(r'\n|\Z', original[start:end]):
                    stop = start + line_end.end()
                    candidates.append(original[:start] + 'writable_roots = ' + value + '\n' + original[stop:])
            else:
                candidates.append(original[:header.end()] + '\nwritable_roots = ' + value + original[header.end():])
        else:
            candidates.append(original + '\n[sandbox_workspace_write]\nwritable_roots = ' + value + '\n')
        for candidate in candidates:
            try:
                if toml.loads(candidate) == expected:
                    updated = candidate
                    break
            except ValueError:
                pass
        else:
            raise ValueError('Cannot safely edit this TOML layout. Add {} to sandbox_workspace_write.writable_roots manually.'.format(state))

    print('Skill: {} -> {}'.format(target, source))
    print('Personal state: {}'.format(state))
    print('Config: {} ({})'.format(config, 'update required' if updated != original else 'already configured'))
    if args.dry_run:
        return
    state.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryFile(dir=str(state)):
        pass
    config.parent.mkdir(parents=True, exist_ok=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    if updated != original:
        if config.exists():
            fd, backup = tempfile.mkstemp(prefix='config.toml.backup-', dir=str(config.parent))
            os.close(fd)
            shutil.copy2(str(config), backup)
            print('Backup: {}'.format(backup))
        fd, temporary = tempfile.mkstemp(prefix='.config.toml-', dir=str(config.parent))
        try:
            with os.fdopen(fd, 'w', encoding='utf-8') as stream:
                stream.write(updated)
            if config.exists():
                shutil.copymode(str(config), temporary)
            os.replace(temporary, str(config))
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    if not target.is_symlink():
        target.symlink_to(source, target_is_directory=True)
    print('Registered. Start a new Codex session and check /skills.')
    print('Managed sandbox policies may also need this state path added by the host.')
    if os.environ.get('PY4SIESTA_SKILL_STATE'):
        print('Keep PY4SIESTA_SKILL_STATE set to this path when starting Codex.')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        print('Registration failed: {}'.format(error), file=sys.stderr)
        sys.exit(1)
