#!/usr/bin/env python3
"""Extract ana.565 tower geometry and render the geometric acceptance study."""
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if len(sys.argv) > 2:
        raise SystemExit('Usage: run.sh [new-production-directory]')
    release = Path(os.environ.get('OFFLINE_MAIN', ''))
    calibration = Path(os.environ.get('CALIBRATIONROOT', ''))
    if release.name != 'ana.565' or not release.is_dir():
        raise SystemExit('Set up ana.565 first; OFFLINE_MAIN must end in /ana.565')
    if calibration.resolve() != (release / 'share/calibrations').resolve():
        raise SystemExit('CALIBRATIONROOT must be the ana.565 share/calibrations directory')
    production = (Path(sys.argv[1]) if len(sys.argv) == 2 else PACKAGE / 'output/calorimeter_acceptance/productions/ana565/tower_surface').resolve()
    if production.exists() and any(production.iterdir()):
        raise SystemExit(f'Refusing to overwrite an existing production: {production}\nSupply a new directory to rerun.')
    results = production / 'results'
    results.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['ROOT_INCLUDE_PATH'] = f'{release}/rootmacros:{release}/include:' + env.get('ROOT_INCLUDE_PATH', '')
    env['LD_LIBRARY_PATH'] = f'{release}/lib:{release}/lib64:' + env.get('LD_LIBRARY_PATH', '')
    sources = list((PACKAGE / 'analysis/calorimeter_acceptance').glob('*.C')) + list(HERE.glob('*.C')) + [HERE / 'run.py', HERE / 'run.sh']
    validator = PACKAGE / 'validation/calorimeter_acceptance/check_outputs.py'
    sources.append(validator)
    sources += [release / 'rootmacros' / name for name in ['G4_CEmc_Spacal.C', 'G4_HcalIn_ref.C', 'G4_HcalOut_ref.C', 'GlobalVariables.C']]
    sources += sorted((calibration / 'CEMC/Geometry_2023ProjTilted').glob('*.xml'))
    for name in ['libg4calo.so', 'libg4detectors.so', 'libg4ihcal.so', 'libg4ohcal.so', 'libphparameter.so', 'libcalo_io.so']:
        candidates = [release / folder / name for folder in ['lib', 'lib64']]
        sources.append(next(path for path in candidates if path.is_file()))
    manifest = {'status': 'running', 'release': str(release), 'calibration_root': str(calibration),
                'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'definition': 'Outer envelope of tower eta-bin boundaries on common tower-center radius',
                'assumptions': {'vertex_x_y_cm': [0, 0], 'full_phi': True, 'internal_eta_and_phi_gaps_ignored': True,
                                'particle_transport': False, 'CDB_or_GDML_loading': False},
                'geometry_method': 'Installed standard setup macros; CEMC XML + macro overrides -> Spacal cell builder -> RawTowerBuilder; HCAL subsystem parameters -> HcalRawTowerBuilder',
                'scan_cm': [-100, 100, 1], 'display_vertices_cm': list(range(0, 101, 10)),
                'source_sha256': {str(p): digest(p) for p in sources}, 'commands': []}
    manifest_path = production / 'manifest.json'
    def save():
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    save()
    try:
        for macro, log in [('ExtractCalorimeterGeometry.C', 'extract.log'), ('PlotCalorimeterAcceptance.C', 'plot.log')]:
            command = ['root', '-l', '-b', '-q', f'{HERE / macro}({json.dumps(str(results))})']
            manifest['commands'].append(command)
            save()
            with (production / log).open('w') as stream:
                subprocess.run(command, cwd=HERE, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
        spec = importlib.util.spec_from_file_location('check_calo_acceptance', validator)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        report = module.check(results)
        (results / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
        manifest['validation'] = report
        manifest['result_sha256'] = {p.name: digest(p) for p in sorted(results.iterdir()) if p.is_file()}
        manifest['status'] = 'complete'
    except BaseException as error:
        manifest['status'] = 'failed'
        manifest['error'] = str(error)
        raise
    finally:
        save()
    print(f'Validated results: {results}')


if __name__ == '__main__':
    main()
