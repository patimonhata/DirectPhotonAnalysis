#!/usr/bin/env python3
"""Build development software, freeze clean releases, and bind new productions to them."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

MODULE = Path(__file__).resolve().parents[1]
REPO = MODULE.parent
DEFAULT_RELEASE = 'ana.565'
STAGES = {name: f'workflows/photon_candidate_selection/run_{name}.sh' for name in ('map', 'reduce', 'merge')}
SOURCE_PATHS = ['PhotonAnalysisTree', 'Pi0Reconstruction/src', 'TopologyEventDisplay', 'README.md', 'AGENTS.md']


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')
    temporary.replace(path)


def git(*args):
    return subprocess.check_output(['git', '-C', str(REPO), *args], text=True).strip()


def source_state():
    paths = git('ls-files', '--cached', '--others', '--exclude-standard', '--', *SOURCE_PATHS).splitlines()
    value = hashlib.sha256()
    for name in sorted(set(paths)):
        # Input lists are data, not the executable source snapshot.
        if name.startswith('PhotonAnalysisTree/input/'):
            continue
        path = REPO / name
        value.update(name.encode() + b'\0' + (digest(path).encode() if path.is_file() else b'missing') + b'\0')
    return {'git_commit': git('rev-parse', 'HEAD'), 'working_tree_dirty': bool(git('status', '--porcelain')), 'source_sha256': value.hexdigest()}


def compiled_state():
    return {str(p.relative_to(REPO)): digest(p) for directory in (MODULE / 'src', REPO / 'Pi0Reconstruction/src')
            for p in sorted(directory.iterdir()) if p.is_file() and (p.suffix in ('.cc', '.h') or p.name == 'CMakeLists.txt')}


def fixed_release(value):
    if not re.fullmatch(r'ana\.[0-9]+', value):
        raise ValueError('Use a fixed release such as ana.565; rolling ana is not accepted')
    return value


def release_environment(release):
    fixed_release(release)
    env = dict(os.environ)
    env.pop('BASH_ENV', None)
    result = subprocess.run(['bash', '--noprofile', '--norc', '-c',
        'set -eo pipefail; source /opt/sphenix/core/bin/sphenix_setup.sh -n "$1" >/dev/null; env -0',
        'photon-software', release], env=env, check=True, capture_output=True)
    selected = dict(item.decode().split('=', 1) for item in result.stdout.split(b'\0') if b'=' in item)
    if Path(selected.get('OFFLINE_MAIN', '')).resolve().name != release:
        raise ValueError('The setup script did not select the requested fixed release')
    return selected


def runtime_environment(prefix, env):
    env = dict(env)
    env['LD_LIBRARY_PATH'] = ':'.join([str(prefix / 'lib64'), str(prefix / 'dependencies/Pi0Reconstruction/lib'), env.get('LD_LIBRARY_PATH', '')])
    env['ROOT_INCLUDE_PATH'] = ':'.join([str(prefix / 'include'), str(prefix / 'dependencies/Pi0Reconstruction/include/pi0reconstruction'), env.get('ROOT_INCLUDE_PATH', '')])
    return env


def inventory(prefix):
    result = {}
    for path in sorted(prefix.rglob('*')):
        if '__pycache__' in path.parts:
            continue
        if path.name in ('manifest.json', 'build-manifest.json') and path.parent == prefix:
            continue
        if path.is_symlink() and not path.resolve().is_relative_to(prefix.resolve()):
            raise ValueError(f'Artifact symlink escapes software directory: {path}')
        if path.is_file():
            result[str(path.relative_to(prefix))] = digest(path)
    return result


def verify(prefix, frozen=True):
    prefix = Path(prefix).resolve()
    manifest = prefix / ('manifest.json' if frozen else 'build-manifest.json')
    value = read_json(manifest)
    if value.get('format_version') != 1 or value.get('kind') != ('frozen' if frozen else 'development'):
        raise ValueError(f'Unexpected software manifest: {manifest}')
    fixed_release(value['environment']['release'])
    if frozen and value['source']['working_tree_dirty']:
        raise ValueError('A production release must have clean source')
    for relative, expected in value['artifacts'].items():
        path = prefix / relative
        if not path.resolve().is_relative_to(prefix) or not path.is_file() or digest(path) != expected:
            raise ValueError(f'Software artifact is missing or changed: {path}')
    if inventory(prefix) != value['artifacts']:
        raise ValueError(f'Software directory contains unrecorded artifacts: {prefix}')
    return value


def cmake_configuration(path):
    values = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith(('#', '//')) or ':' not in line or '=' not in line:
            continue
        key, value = line.split('=', 1)
        name = key.split(':', 1)[0]
        if name.startswith(('CMAKE_CXX_', 'CMAKE_SHARED_LINKER_FLAGS')) or name in ('CMAKE_BUILD_TYPE', 'ROOT_DIR', 'PHOTON_CONFIGURED_PREFIX'):
            values[name] = value
    return values


def build(args):
    release = fixed_release(args.release)
    env = release_environment(release)
    build_dir = (args.build_dir or MODULE / 'build/dev' / release).resolve()
    prefix = (args.install_prefix or MODULE / 'install/dev' / release).resolve()
    if prefix == (MODULE / 'install').resolve() or prefix.is_relative_to((MODULE / 'software').resolve()):
        raise ValueError('Development builds must not overwrite the legacy install or frozen software')
    if (prefix / 'manifest.json').exists():
        raise ValueError('Refusing to update a frozen software directory')
    prefix.mkdir(parents=True, exist_ok=True)
    # Invalidate the previous record before any build/install mutation.
    (prefix / 'build-manifest.json').unlink(missing_ok=True)
    compiled_before = compiled_state()
    dependency = prefix / 'dependencies/Pi0Reconstruction'
    compiler = shutil.which('g++', path=env['PATH'])
    common = ['-DCMAKE_BUILD_TYPE=RelWithDebInfo', f'-DCMAKE_CXX_COMPILER={compiler}']
    configs = [
        (REPO / 'Pi0Reconstruction/src', build_dir / 'pi0reconstruction', dependency, []),
        (MODULE / 'src', build_dir / 'photon', prefix, [f'-DPI0RECONSTRUCTION_INCLUDE_DIR={dependency}/include/pi0reconstruction',
             f'-DPI0RECONSTRUCTION_LIBRARY={dependency}/lib/libPi0Reconstruction.so']),
    ]
    for source, directory, destination, extra in configs:
        cache = directory / 'CMakeCache.txt'
        if cache.exists():
            values = {line.split(':', 1)[0]: line.split('=', 1)[1] for line in cache.read_text().splitlines() if ':' in line and '=' in line and not line.startswith('//')}
            if values.get('CMAKE_HOME_DIRECTORY') != str(source) or values.get('CMAKE_INSTALL_PREFIX') != str(destination):
                raise ValueError(f'Build cache belongs to another source/install; use a new build directory: {directory}')
            previous = values.get('PHOTON_CONFIGURED_PREFIX')
            if previous and previous != env['OFFLINE_MAIN']:
                raise ValueError(f'Build cache belongs to another release: {directory}')
        subprocess.run(['cmake', '-S', str(source), '-B', str(directory), f'-DCMAKE_INSTALL_PREFIX={destination}', *common, *extra], env=env, check=True)
        subprocess.run(['cmake', '--build', str(directory), '--parallel', str(args.jobs)], env=env, check=True)
        subprocess.run(['cmake', '--install', str(directory)], env=env, check=True)
    if compiled_state() != compiled_before:
        raise ValueError('Compiled source changed during the build; rebuild before recording provenance')
    before = source_state()
    active = runtime_environment(prefix, env)
    resolved = {}
    for name in ('libPhotonAnalysisTree.so', 'libPi0GammaOnnx.so'):
        result = subprocess.run(['ldd', str(prefix / 'lib64' / name)], env=active, check=True, capture_output=True, text=True)
        if 'not found' in result.stdout:
            raise ValueError(result.stdout)
        paths = re.findall(r'=> (/\S+)', result.stdout)
        for path in paths:
            real = Path(path).resolve()
            if real.is_relative_to(REPO) and not real.is_relative_to(prefix):
                raise ValueError(f'Build uses an unfrozen local dependency: {real}')
        resolved[name] = paths
    value = {'format_version': 1, 'kind': 'development', 'source': before,
        'environment': {'release': release, 'offline_main': env['OFFLINE_MAIN'],
            'root_version': subprocess.check_output(['root-config', '--version'], env=env, text=True).strip(),
            'compiler': compiler, 'compiler_version': subprocess.check_output([compiler, '--version'], env=env, text=True).splitlines()[0]},
        'compiled_sources': compiled_before,
        'build': {'directory': str(build_dir), 'install_prefix': str(prefix), 'type': 'RelWithDebInfo', 'resolved_dependencies': resolved,
            'cmake': {directory.name: cmake_configuration(directory / 'CMakeCache.txt') for _, directory, _, _ in configs}},
        'built_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'artifacts': inventory(prefix)}
    write_json(prefix / 'build-manifest.json', value)
    print(f'Development install: {prefix}\nManifest: {prefix / "build-manifest.json"}')


def freeze(args):
    prefix = (args.install_prefix or MODULE / 'install/dev' / args.release).resolve()
    value = verify(prefix, frozen=False)
    if value['environment']['release'] != args.release:
        raise ValueError('Requested release differs from the build being frozen')
    current = source_state()
    if current['working_tree_dirty'] or value['source']['working_tree_dirty'] or current != value['source']:
        raise ValueError('Freeze requires clean source matching the build record; commit changes and rebuild first')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.software_id):
        raise ValueError('Invalid software ID')
    destination = (args.software_directory or MODULE / 'software') / args.software_id
    if destination.exists() or destination.is_symlink():
        raise ValueError(f'Software ID already exists: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.parent / ('.' + args.software_id + '.lock')
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.close(descriptor)
    temporary = None
    try:
        temporary = Path(tempfile.mkdtemp(prefix='.freeze-', dir=destination.parent))
        shutil.copytree(prefix, temporary, dirs_exist_ok=True, symlinks=True)
        snapshot = temporary / 'source'
        snapshot.mkdir()
        archive = temporary / 'source.tar'
        with archive.open('wb') as stream:
            subprocess.run(['git', '-C', str(REPO), 'archive', 'HEAD', *SOURCE_PATHS], check=True, stdout=stream)
        with tarfile.open(archive) as stream:
            # Do not duplicate input lists; map arguments record the external manifest.
            members = [m for m in stream.getmembers() if not m.name.startswith('PhotonAnalysisTree/input/')]
            stream.extractall(snapshot, members=members, filter='data')
        archive.unlink()
        (temporary / 'build-manifest.json').unlink()
        value.update(kind='frozen', software_id=args.software_id, features=['external_map_root_v1'], frozen_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        value['artifacts'] = inventory(temporary)
        write_json(temporary / 'manifest.json', value)
        verify(temporary)
        # mkdir above never reuses an existing release; publish the complete verified tree.
        if destination.exists():
            raise ValueError(f'Software ID already exists: {destination}')
        temporary.rename(destination)
    finally:
        if temporary is not None and temporary.exists():
            shutil.rmtree(temporary)
        lock.unlink()
    print(destination.resolve())


def production(args):
    prefix = args.software.resolve()
    software = verify(prefix)
    manifest = args.manifest.resolve()
    if manifest.exists() or (manifest.parent.exists() and any(manifest.parent.iterdir())):
        raise ValueError('Create a production in a new, empty configuration directory; retained results are never relabelled')
    map_root = getattr(args, 'map_root', None)
    if map_root is not None:
        if 'external_map_root_v1' not in software.get('features', []):
            raise ValueError('Freeze new software with external map-root support before using --map-root')
        map_root = map_root.resolve()
        if map_root.is_relative_to(manifest.parent) or manifest.parent.is_relative_to(map_root):
            raise ValueError('--map-root must be separate from the production configuration directory')
        if map_root.exists() and (not map_root.is_dir() or any(map_root.iterdir())):
            raise ValueError('--map-root must be a new or empty directory; retained maps are never relabelled')
    value = {'format_version': 1, 'workflow': 'photon_candidate_selection', 'production': args.production,
        'configuration': args.configuration, 'status': 'configured',
        'software': {'id': software['software_id'], 'path': str(prefix), 'manifest_sha256': digest(prefix / 'manifest.json'),
            'git_commit': software['source']['git_commit'], 'release': software['environment']['release']},
        'settings': read_json(args.settings) if args.settings else {},
        'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    if map_root is not None:
        value['output_roots'] = {'map': str(map_root)}
        map_root.mkdir(parents=True, exist_ok=True)
        manifest.parent.mkdir(parents=True, exist_ok=True)
        link = manifest.parent / 'maps'
        link.symlink_to(map_root, target_is_directory=True)
        try:
            write_json(manifest, value)
        except Exception:
            link.unlink()
            raise
    else:
        write_json(manifest, value)
    print(manifest)


def selection(software=None, production_manifest=None, release=None):
    production_path = Path(production_manifest).resolve() if production_manifest else None
    if production_path:
        value = read_json(production_path)
        if value.get('format_version') != 1 or value.get('workflow') != 'photon_candidate_selection':
            raise ValueError('This production has no supported software binding')
        bound = value['software']
        if software and Path(software).resolve() != Path(bound['path']).resolve():
            raise ValueError('Software override differs from the production binding')
        software = bound['path']
        if digest(Path(software) / 'manifest.json') != bound['manifest_sha256']:
            raise ValueError('Production software manifest has changed')
    if software:
        prefix = Path(software).resolve()
        value = verify(prefix)
        if production_path and (value['software_id'] != bound['id'] or value['source']['git_commit'] != bound['git_commit'] or value['environment']['release'] != bound['release']):
            raise ValueError('Production software identity does not match its manifest')
        source = prefix / 'source/PhotonAnalysisTree'
    else:
        selected = fixed_release(release or DEFAULT_RELEASE)
        prefix = Path(os.environ.get('PHOTON_INSTALL_PREFIX', MODULE / 'install/dev' / selected)).resolve()
        value = verify(prefix, frozen=False)
        source = MODULE
    if release and value['environment']['release'] != release:
        raise ValueError('Requested release differs from the software build')
    return prefix, source, value, production_path


def resolve(args):
    prefix, source, value, _ = selection(args.software, args.production_manifest, args.release)
    for item in (value['environment']['release'], prefix, source):
        if '\n' in str(item):
            raise ValueError('Runtime paths must not contain newlines')
        print(item)


def assert_output(path, production_path, stage=None):
    if not production_path:
        return
    production_path = Path(production_path).resolve()
    resolved = Path(path).resolve()
    if resolved.is_relative_to(production_path.parent):
        return
    # Only map jobs may write to the explicitly registered external root. Keep the
    # stored canonical root fixed so retargeting a symlink cannot expand this grant.
    root = read_json(production_path).get('output_roots', {}).get('map') if stage == 'map' else None
    if root and Path(root).is_absolute() and resolved.is_relative_to(Path(root)):
        return
    raise ValueError('Output lies outside the bound production configuration or registered map root for this stage')


def check_inputs(args):
    import glob
    prefix, _, value, production_path = selection(args.software, args.production_manifest, args.release)
    if not production_path:
        return
    matches = sorted(glob.glob(args.pattern))
    if not matches:
        raise ValueError('No inputs match the production stage')
    for name in matches:
        path = Path(name).resolve()
        assert_output(path, production_path, args.stage)
        record = read_json(str(path) + '.execution.json')
        if (record.get('stage') != args.stage or record.get('software_manifest_sha256') != digest(prefix / 'manifest.json')
                or record.get('production_manifest_sha256') != digest(production_path) or record.get('output_sha256') != digest(path)):
            raise ValueError(f'Input has different or invalid production provenance: {path}')


def record(args):
    prefix, _, value, production_path = selection(args.software, args.production_manifest, args.release)
    output = args.output.resolve()
    assert_output(output, production_path, args.stage)
    if args.final_output is not None:
        assert_output(args.final_output, production_path, args.stage)
    manifest = prefix / ('manifest.json' if value['kind'] == 'frozen' else 'build-manifest.json')
    selected_release = os.environ.get('OFFLINE_MAIN', '')
    if Path(selected_release).resolve() != Path(value['environment']['offline_main']).resolve():
        raise ValueError('Runtime release differs from the recorded build environment')
    provenance = {'format_version': 1, 'stage': args.stage, 'arguments': args.arguments, 'output': str((args.final_output or output).resolve()),
        'output_sha256': digest(output), 'software_id': value.get('software_id', 'development'),
        'software_manifest_sha256': digest(manifest), 'source': value['source'],
        'runtime_source': source_state() if value['kind'] == 'development' else value['source'],
        'environment': value['environment'], 'runtime_offline_main': selected_release,
        'production_manifest': str(production_path) if production_path else None,
        'production_manifest_sha256': digest(production_path) if production_path else None}
    provenance['completed_at_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    if args.stage == 'map' and len(args.arguments) >= 5:
        input_manifest = Path(args.arguments[4]).resolve()
        provenance['input_manifest'] = {'path': str(input_manifest), 'sha256': digest(input_manifest)}
        expected = os.environ.get('PHOTON_INPUT_MANIFEST_SHA256')
        if expected and expected != provenance['input_manifest']['sha256']:
            raise ValueError('Input manifest changed during map generation')
    sidecar = Path(str(output) + '.execution.json')
    if sidecar.exists():
        raise ValueError(f'Execution record already exists: {sidecar}')
    write_json(sidecar, provenance)


def run(args):
    prefix, source, value, production_path = selection(args.software, args.production_manifest, args.release)
    env = dict(os.environ)
    env['PHOTON_SOFTWARE_ROOT'] = str(prefix) if value['kind'] == 'frozen' else ''
    if production_path:
        env['PHOTON_PRODUCTION_MANIFEST'] = str(production_path)
    env.pop('BASH_ENV', None)
    os.execve(source / STAGES[args.stage], [str(source / STAGES[args.stage]), *args.arguments], env)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build', help='build both local libraries in a fixed release')
    b.add_argument('--release', type=fixed_release, default=DEFAULT_RELEASE)
    b.add_argument('--build-dir', type=Path)
    b.add_argument('--install-prefix', type=Path)
    b.add_argument('--jobs', type=int, default=int(os.environ.get('BUILD_JOBS', '4')))
    f = sub.add_parser('freeze', help='publish a new software ID from a clean matching build')
    f.add_argument('software_id')
    f.add_argument('--release', type=fixed_release, default=DEFAULT_RELEASE)
    f.add_argument('--install-prefix', type=Path)
    f.add_argument('--software-directory', type=Path)
    v = sub.add_parser('verify', help='verify every frozen artifact and source file')
    v.add_argument('software', type=Path)
    init = sub.add_parser('production', help='bind a new empty production configuration to frozen software')
    init.add_argument('software', type=Path)
    init.add_argument('manifest', type=Path)
    init.add_argument('--production', required=True)
    init.add_argument('--configuration', required=True)
    init.add_argument('--settings', type=Path)
    init.add_argument('--map-root', type=Path, help='new or empty external map directory; creates configuration/maps as a symlink')
    for name in ('resolve', 'check-inputs', 'record', 'run', 'check-output'):
        command = sub.add_parser(name)
        command.add_argument('--software', default=os.environ.get('PHOTON_SOFTWARE_ROOT') or None)
        command.add_argument('--production-manifest', default=os.environ.get('PHOTON_PRODUCTION_MANIFEST') or None)
        command.add_argument('--release', type=fixed_release, default=os.environ.get('PHOTON_RELEASE') or None)
        if name == 'check-inputs':
            command.add_argument('stage', choices=STAGES)
            command.add_argument('pattern')
        if name in ('record', 'check-output'):
            if name == 'record':
                command.add_argument('--final-output', type=Path)
                command.add_argument('stage', choices=STAGES)
            if name == 'check-output':
                command.add_argument('--stage', choices=STAGES)
            command.add_argument('output', type=Path)
        if name == 'run':
            command.add_argument('stage', choices=STAGES)
        if name in ('record', 'run'):
            command.add_argument('arguments', nargs=argparse.REMAINDER)
    return p


def main():
    args = parser().parse_args()
    try:
        if args.command == 'verify':
            value = verify(args.software)
            print(f'Verified {value["software_id"]}: {len(value["artifacts"])} artifacts')
        elif args.command == 'check-output':
            _, _, _, production_path = selection(args.software, args.production_manifest, args.release)
            assert_output(args.output, production_path, args.stage)
        else:
            globals()[args.command.replace('-', '_')](args)
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
        print(f'photon software: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
