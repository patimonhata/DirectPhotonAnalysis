#!/usr/bin/env python3
"""Validate geometry snapshots and PDFs using independent angular calculations."""
import csv
import json
import math
from pathlib import Path
import re
import sys


def check(results):
    results = Path(results)
    def read(name):
        with (results / name).open() as stream:
            return list(csv.DictReader(stream))
    geometry = {r['detector']: r for r in read('geometry.csv')}
    assert set(geometry) == {'CEMC', 'HCALIN', 'HCALOUT'}
    bins = read('eta_bins.csv')
    gap_report = {}
    for name, geo in geometry.items():
        rows = [r for r in bins if r['detector'] == name]
        assert len(rows) == int(geo['eta_bins'])
        assert [int(r['bin']) for r in rows] == list(range(len(rows)))
        assert math.isclose(min(float(r['eta_low']) for r in rows), float(geo['eta_min']), abs_tol=1e-12)
        assert math.isclose(max(float(r['eta_high']) for r in rows), float(geo['eta_max']), abs_tol=1e-12)
        gap_report[name] = max([0.] + [float(b['eta_low']) - float(a['eta_high']) for a, b in zip(rows, rows[1:])])
    checked = 0
    for file, vertices in [('acceptance_scan.csv', range(-100, 101)), ('acceptance_selected_vertices.csv', range(0, 101, 10))]:
        rows = read(file)
        assert len(rows) == 3 * len(vertices)
        for name, geo in geometry.items():
            subset = [r for r in rows if r['detector'] == name]
            assert [int(r['z_vertex_cm']) for r in subset] == list(vertices)
            radius = float(geo['radius_cm'])
            previous = None
            lookup = {int(r['z_vertex_cm']): r for r in subset}
            for row in subset:
                vertex = int(row['z_vertex_cm'])
                lo, hi = float(row['eta_min']), float(row['eta_max'])
                assert math.isfinite(lo) and math.isfinite(hi) and lo < hi
                for key, eta in [('z_min_cm', lo), ('z_max_cm', hi)]:
                    endpoint = float(geo[key])
                    theta = math.atan2(radius, endpoint - vertex)
                    assert abs(eta + math.log(math.tan(theta / 2))) < 1e-12
                    assert abs(vertex + radius * math.sinh(eta) - endpoint) < 1e-9
                assert abs(float(row['delta_eta']) - (hi - lo)) < 1e-12
                if previous:
                    assert lo < previous[0] and hi < previous[1]
                previous = (lo, hi)
                if -vertex in lookup:
                    assert abs(lo + float(lookup[-vertex]['eta_max'])) < 1e-12
                checked += 1
    pages = {}
    for name, expected in [('acceptance_vs_z.pdf', 1), ('detector_comparison.pdf', 1), ('acceptance_vertex_views.pdf', 11)]:
        data = (results / name).read_bytes()
        assert data.startswith(b'%PDF-') and b'%%EOF' in data[-1024:]
        pages[name] = len(re.findall(rb'/Type\s*/Page\b', data))
        assert pages[name] == expected, (name, pages[name])
    assert (results / 'tower_geometry.root').read_bytes()[:4] == b'root'
    return {'status': 'passed', 'checked_rows': checked, 'pdf_pages': pages,
            'largest_ignored_internal_eta_gap': gap_report,
            'checks': ['theta-based eta', 'ray endpoint recovery', 'vertex monotonicity', 'reflection symmetry',
                       'bin edge envelope', 'complete scan grids', 'PDF structure and page counts']}


if __name__ == '__main__':
    report = check(sys.argv[1])
    print(json.dumps(report, indent=2))
