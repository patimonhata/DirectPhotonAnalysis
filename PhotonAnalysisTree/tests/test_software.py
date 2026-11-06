"""Regression tests for frozen artifacts, clean-source gates, and production lineage."""
import argparse
import importlib.util
import json
import os
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

module = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('photon_software', module / 'tools/software.py')
software = importlib.util.module_from_spec(spec)
spec.loader.exec_module(software)


class SoftwareTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='photon-software-test-')
        self.root = Path(self.temporary.name)
        self.prefix = self.root / 'software/test'
        (self.prefix / 'lib64').mkdir(parents=True)
        (self.prefix / 'lib64/libPhotonAnalysisTree.so').write_bytes(b'library v1')
        self.source = {'git_commit': 'abc123', 'working_tree_dirty': False, 'source_sha256': 'source v1'}
        self.manifest = {'format_version': 1, 'kind': 'frozen', 'software_id': 'test', 'source': self.source,
                         'environment': {'release': 'ana.565', 'offline_main': '/release/ana.565'},
                         'artifacts': software.inventory(self.prefix)}
        software.write_json(self.prefix / 'manifest.json', self.manifest)
        self.production = self.root / 'output/new/configuration/manifest.json'
        software.production(argparse.Namespace(software=self.prefix, manifest=self.production, production='new', configuration='configuration', settings=None))

    def tearDown(self):
        self.temporary.cleanup()

    def enable_external_maps(self):
        self.manifest['features'] = ['external_map_root_v1']
        software.write_json(self.prefix / 'manifest.json', self.manifest)
        self.production = self.root / 'user/external/configuration/manifest.json'
        self.map_root = self.root / 'tg/external/maps'
        return argparse.Namespace(software=self.prefix, manifest=self.production, production='external',
                                  configuration='configuration', settings=None, map_root=self.map_root)

    def test_external_map_root_creation_and_stage_limits(self):
        software.production(self.enable_external_maps())
        local = self.production.parent
        self.assertTrue((local / 'maps').is_symlink())
        self.assertEqual((local / 'maps').resolve(), self.map_root.resolve())
        self.assertEqual(software.read_json(self.production)['output_roots'], {'map': str(self.map_root.resolve())})
        for path in (local / 'maps/jet3/map_000000.root', self.map_root / 'jet3/.temporary.root'):
            software.assert_output(path, self.production, 'map')
            for stage in (None, 'reduce', 'merge'):
                with self.assertRaisesRegex(ValueError, 'outside'):
                    software.assert_output(path, self.production, stage)
        software.assert_output(local / 'partials/jet/jet3/shard_0/result.root', self.production, 'reduce')
        software.assert_output(local / 'results/jet/result.root', self.production, 'merge')
        for path in (self.root / 'tg/unregistered/map.root', self.map_root.with_name('maps_suffix') / 'map.root'):
            with self.assertRaisesRegex(ValueError, 'outside'):
                software.assert_output(path, self.production, 'map')
        escape = self.root / 'tg/other'
        escape.mkdir()
        (self.map_root / 'escape').symlink_to(escape, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'outside'):
            software.assert_output(local / 'maps/escape/map.root', self.production, 'map')
        (local / 'maps').unlink()
        (local / 'maps').symlink_to(escape, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'outside'):
            software.assert_output(local / 'maps/map.root', self.production, 'map')
        self.map_root.rename(self.map_root.with_name('original_maps'))
        self.map_root.symlink_to(escape, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'outside'):
            software.assert_output(self.map_root / 'map.root', self.production, 'map')

    def test_external_map_root_requires_new_software_and_empty_separate_directory(self):
        args = argparse.Namespace(software=self.prefix, manifest=self.root / 'new/manifest.json', production='new',
                                  configuration='configuration', settings=None, map_root=self.root / 'tg/maps')
        with self.assertRaisesRegex(ValueError, 'Freeze new software'):
            software.production(args)
        args = self.enable_external_maps()
        self.map_root.mkdir(parents=True)
        existing = self.map_root / 'retained.root'
        existing.write_bytes(b'retained')
        with self.assertRaisesRegex(ValueError, 'empty directory'):
            software.production(args)
        self.assertEqual(existing.read_bytes(), b'retained')
        self.assertFalse(self.production.exists())
        for path in (self.production.parent, self.production.parent / 'maps', self.production.parent.parent):
            args.map_root = path
            with self.assertRaisesRegex(ValueError, 'separate'):
                software.production(args)
        args.map_root = existing
        with self.assertRaisesRegex(ValueError, 'empty directory'):
            software.production(args)
        existing.unlink()
        args.map_root = self.map_root
        software.production(args)  # Pre-created empty storage directories are accepted.

    def test_external_map_roundtrip_and_provenance_rejection(self):
        software.production(self.enable_external_maps())
        temporary = self.production.parent / 'maps/jet3/.temporary.root'
        temporary.parent.mkdir()
        temporary.write_bytes(b'map contents')
        final = temporary.with_name('map_000000.root')
        args = argparse.Namespace(software=self.prefix, production_manifest=self.production, release=None,
                                  stage='map', output=temporary, final_output=final, arguments=[])
        with patch.dict(os.environ, {'OFFLINE_MAIN': '/release/ana.565'}):
            software.record(args)
        temporary.rename(final)
        Path(str(temporary) + '.execution.json').rename(str(final) + '.execution.json')
        saved = software.read_json(str(final) + '.execution.json')
        self.assertEqual(saved['output'], str(final.resolve()))
        self.assertEqual(saved['production_manifest'], str(self.production.resolve()))
        inputs = argparse.Namespace(software=self.prefix, production_manifest=self.production, release=None,
                                    stage='map', pattern=str(final.parent / 'map_*.root'))
        software.check_inputs(inputs)
        inputs.pattern = str(final.resolve())
        software.check_inputs(inputs)
        # A local reducer can consume external maps and retain its own output locally.
        partial = self.production.parent / 'partials/jet/jet3/shard_0/partial.root'
        partial.parent.mkdir(parents=True)
        partial.write_bytes(b'partial contents')
        args.stage = 'reduce'; args.output = partial; args.final_output = None
        with patch.dict(os.environ, {'OFFLINE_MAIN': '/release/ana.565'}):
            software.record(args)
        inputs.stage = 'reduce'; inputs.pattern = str(partial)
        software.check_inputs(inputs)
        inputs.stage = 'map'; inputs.pattern = str(final)
        saved['production_manifest_sha256'] = 'another production'
        software.write_json(str(final) + '.execution.json', saved)
        with self.assertRaisesRegex(ValueError, 'provenance'):
            software.check_inputs(inputs)
        # Valid-looking hashes cannot authorize an unregistered location.
        saved['production_manifest_sha256'] = software.digest(self.production)
        outside = self.root / 'tg/unregistered.root'
        outside.write_bytes(final.read_bytes())
        software.write_json(str(outside) + '.execution.json', saved)
        inputs.pattern = str(outside)
        with self.assertRaisesRegex(ValueError, 'outside'):
            software.check_inputs(inputs)
        args.stage = 'map'; args.output = final; args.final_output = outside
        with patch.dict(os.environ, {'OFFLINE_MAIN': '/release/ana.565'}):
            with self.assertRaisesRegex(ValueError, 'outside'):
                software.record(args)

    def test_external_map_check_output_shell_forwards_stage(self):
        software.production(self.enable_external_maps())
        env = {k: v for k, v in os.environ.items() if not k.startswith('PHOTON_') and k != 'BASH_ENV'}
        env.update(PHOTON_SOURCE_ROOT=str(module), PHOTON_PRODUCTION_MANIFEST=str(self.production))
        command = ['bash', '--noprofile', '--norc', '-c',
                   'set -euo pipefail; source "$1"; photon_check_output "$2" --stage "$3"', 'output-test',
                   str(module / 'workflows/common/runtime.sh'), str(self.production.parent / 'maps/jet3/map.root')]
        for stage, expected in (('map', 0), ('reduce', 1), ('merge', 1)):
            result = subprocess.run([*command, stage], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, expected, result.stderr)

    def test_rolling_release_is_rejected(self):
        self.assertEqual(software.fixed_release('ana.565'), 'ana.565')
        for value in ('ana', 'ana.latest', '../ana.565'):
            with self.assertRaises(ValueError):
                software.fixed_release(value)

    def test_artifact_tampering_and_extra_files_are_rejected(self):
        software.verify(self.prefix)
        library = self.prefix / 'lib64/libPhotonAnalysisTree.so'
        library.write_bytes(b'library v2')
        with self.assertRaisesRegex(ValueError, 'changed'):
            software.verify(self.prefix)
        library.write_bytes(b'library v1')
        (self.prefix / 'lib64/extra.so').write_bytes(b'extra')
        with self.assertRaisesRegex(ValueError, 'unrecorded'):
            software.verify(self.prefix)

    def test_artifact_cannot_escape_bundle(self):
        self.manifest['artifacts'] = {'../outside': 'anything'}
        software.write_json(self.prefix / 'manifest.json', self.manifest)
        with self.assertRaises(ValueError):
            software.verify(self.prefix)

    def test_existing_results_are_not_relabelled(self):
        args = argparse.Namespace(software=self.prefix, manifest=self.production, production='other', configuration='other', settings=None)
        with self.assertRaisesRegex(ValueError, 'empty'):
            software.production(args)
        self.production.unlink()
        (self.production.parent / 'retained.root').write_bytes(b'retained')
        with self.assertRaisesRegex(ValueError, 'empty'):
            software.production(args)

    def test_production_rejects_changed_manifest_or_override(self):
        software.selection(production_manifest=self.production)
        with self.assertRaisesRegex(ValueError, 'override'):
            software.selection(software=self.root / 'other', production_manifest=self.production)
        self.manifest['source']['git_commit'] = 'different'
        software.write_json(self.prefix / 'manifest.json', self.manifest)
        with self.assertRaisesRegex(ValueError, 'changed'):
            software.selection(production_manifest=self.production)

    def test_runtime_dispatch_uses_frozen_entrypoint(self):
        entry = self.prefix / 'source/PhotonAnalysisTree/workflows/photon_candidate_selection/run_map.sh'
        entry.parent.mkdir(parents=True)
        entry.write_text('printf "frozen entry: %s\\n" "$1"\n')
        self.manifest['artifacts'] = software.inventory(self.prefix)
        software.write_json(self.prefix / 'manifest.json', self.manifest)
        self.production.unlink()
        software.production(argparse.Namespace(software=self.prefix, manifest=self.production, production='new', configuration='configuration', settings=None))
        env = {k: v for k, v in os.environ.items() if not k.startswith('PHOTON_') and k != 'BASH_ENV'}
        env['PHOTON_PRODUCTION_MANIFEST'] = str(self.production)
        result = subprocess.run(['bash', '--noprofile', '--norc', '-c',
            'set -euo pipefail; source "$1"; photon_dispatch "$2" workflows/photon_candidate_selection/run_map.sh example; exit 91',
            'dispatch-test', str(module / 'workflows/common/runtime.sh'), str(module)], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), 'frozen entry: example')

    def test_freeze_requires_matching_clean_build(self):
        args = argparse.Namespace(install_prefix=self.prefix, release='ana.565', software_id='another', software_directory=self.root / 'software')
        with patch.object(software, 'verify', return_value=self.manifest):
            for source in ({**self.source, 'working_tree_dirty': True}, {**self.source, 'git_commit': 'different'}):
                with patch.object(software, 'source_state', return_value=source):
                    with self.assertRaisesRegex(ValueError, 'clean source'):
                        software.freeze(args)
        self.assertFalse((self.root / 'software/another').exists())

    def test_software_id_is_not_overwritten(self):
        args = argparse.Namespace(install_prefix=self.prefix, release='ana.565', software_id='test', software_directory=self.root / 'software')
        original = software.digest(self.prefix / 'manifest.json')
        with patch.object(software, 'verify', return_value=self.manifest), patch.object(software, 'source_state', return_value=self.source):
            with self.assertRaisesRegex(ValueError, 'already exists'):
                software.freeze(args)
        self.assertEqual(software.digest(self.prefix / 'manifest.json'), original)

    def test_record_tracks_final_map_and_rejects_runtime_mismatch(self):
        output = self.production.parent / 'maps/jet5/.temporary.root'
        output.parent.mkdir(parents=True)
        output.write_bytes(b'map contents')
        final = output.with_name('map_000000.root')
        args = argparse.Namespace(software=self.prefix, production_manifest=self.production, release=None,
                                  stage='map', output=output, final_output=final, arguments=['jet5', 'cuts'])
        with patch.dict(os.environ, {'OFFLINE_MAIN': '/release/ana.579'}):
            with self.assertRaisesRegex(ValueError, 'Runtime release'):
                software.record(args)
        with patch.dict(os.environ, {'OFFLINE_MAIN': '/release/ana.565'}):
            software.record(args)
        record = software.read_json(str(output) + '.execution.json')
        self.assertEqual(record['output'], str(final.resolve()))
        self.assertEqual(record['output_sha256'], software.digest(output))
        self.assertEqual(record['arguments'], ['jet5', 'cuts'])
        self.assertEqual(record['production_manifest_sha256'], software.digest(self.production))

    def test_stage_inputs_reject_missing_mixed_and_changed_outputs(self):
        output = self.production.parent / 'maps/jet5/map_000000.root'
        output.parent.mkdir(parents=True)
        output.write_bytes(b'map contents')
        args = argparse.Namespace(software=self.prefix, production_manifest=self.production, release=None,
                                  stage='map', pattern=str(output))
        with self.assertRaises(OSError):
            software.check_inputs(args)
        record = {'stage': 'map', 'software_manifest_sha256': software.digest(self.prefix / 'manifest.json'),
                  'production_manifest_sha256': software.digest(self.production), 'output_sha256': software.digest(output)}
        software.write_json(str(output) + '.execution.json', record)
        software.check_inputs(args)
        record['software_manifest_sha256'] = 'another software'
        software.write_json(str(output) + '.execution.json', record)
        with self.assertRaisesRegex(ValueError, 'provenance'):
            software.check_inputs(args)
        record['software_manifest_sha256'] = software.digest(self.prefix / 'manifest.json')
        software.write_json(str(output) + '.execution.json', record)
        output.write_bytes(b'changed after recording')
        with self.assertRaisesRegex(ValueError, 'provenance'):
            software.check_inputs(args)
        with self.assertRaisesRegex(ValueError, 'outside'):
            software.assert_output(self.root / 'retained.root', self.production)


if __name__ == '__main__':
    unittest.main()
