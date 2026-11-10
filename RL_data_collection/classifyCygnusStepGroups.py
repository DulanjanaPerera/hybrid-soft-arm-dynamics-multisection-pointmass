"""Group saved tracker step responses by long rests; no device/pressure I/O."""
from pathlib import Path
import argparse
import csv
import json
import math
import statistics


def xyz(row, prefix='tip'):
    return [float(row[prefix+'_'+n]) for n in ['x_mm', 'y_mm', 'z_mm']]


def baseline_quaternion(rows):
    samples = [[float(r['base_'+n]) for n in ['qw', 'qx', 'qy', 'qz']]
               for r in rows if r['base_status'] == '1' and 1 <= float(r['elapsed_s']) <= 5]
    aligned = [[(-v if sum(a*b for a, b in zip(q, samples[0])) < 0 else v) for v in q] for q in samples]
    q = [statistics.mean(row[j] for row in aligned) for j in range(4)]
    norm = math.sqrt(sum(v*v for v in q))
    return [v/norm for v in q]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('folder', type=Path)
    ap.add_argument('--order', nargs='+', default=['P3', 'P2', 'P1'])
    ap.add_argument('--compare', type=Path)
    args = ap.parse_args()
    rows = list(csv.DictReader((args.folder/'tracking.csv').open()))
    times = [float(r['elapsed_s']) for r in rows]
    baseline = [statistics.median(xyz(r)[j] for r in rows if r['tip_status'] == '1' and 1 <= float(r['elapsed_s']) <= 5) for j in range(3)]
    active = [i for i, r in enumerate(rows) if r['tip_status'] == '1' and r['base_status'] == '1' and math.dist(xyz(r), baseline) > 30]
    groups = []
    for i in active:
        if not groups or times[i]-times[groups[-1][-1]] > 8:
            groups.append([])
        groups[-1].append(i)
    assert len(groups) == len(args.order), (len(groups), args.order)
    boundaries = [0]+[(times[a[-1]]+times[b[0]])/2 for a, b in zip(groups, groups[1:])]+[times[-1]+1e-6]
    records = []
    for j, group in enumerate(groups):
        subset = [r for r in rows if boundaries[j] <= float(r['elapsed_s']) < boundaries[j+1] and r['startup'] == '0']
        rec = dict(PMA=args.order[j], response_start_s=times[group[0]], response_end_s=times[group[-1]],
                   group_start_s=boundaries[j], group_end_s=min(boundaries[j+1], times[-1]), runtime_reports=len(subset),
                   tip_missing_reports=sum(r['tip_status'] == '2' for r in subset),
                   base_missing_reports=sum(r['base_status'] == '2' for r in subset))
        for axis, k in zip('XYZ', range(3)):
            rec['median_raw_tip_deflection_'+axis+'_mm'] = statistics.median(xyz(rows[i])[k] for i in group)-baseline[k]
        records.append(rec)
    with (args.folder/'movement_groups.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    result = dict(method='Valid raw-tip displacement >30mm from1-5sbaseline; join active samples across <=8s; split long rest gaps at midpoints.',
                  PMA_order_source='User-confirmed activation order; no pressure timestamps in tracker log', groups=records)
    if args.compare:
        old = list(csv.DictReader((args.compare/'tracking.csv').open()))
        qa, qb = baseline_quaternion(old), baseline_quaternion(rows)
        rotation = 2*math.acos(min(1, abs(sum(a*b for a, b in zip(qa, qb)))))*180/math.pi
        result['comparison'] = dict(previous=str(args.compare), baseline_base_sensor_orientation_change_deg=rotation)
        for name, data in [('previous', old), ('current', rows)]:
            indices = sorted({j for i, r in enumerate(data) if r['startup'] == '0' and r['tip_status'] == '2'
                              for j in [i-1, i+1] if 0 <= j < len(data) and data[j]['tip_status'] == '1' and data[j]['base_status'] == '1'})
            result['comparison'][name+'_adjacent_valid_tip_Z_range_mm'] = ([min(xyz(data[i])[2] for i in indices),
                                                                          statistics.median(xyz(data[i])[2] for i in indices),
                                                                          max(xyz(data[i])[2] for i in indices)] if indices else [])
        result['comparison']['caveat'] = 'Adjacent valid coordinates are not measured coordinates during a missing transform; raw tracker-frame positions, no sensor-offset correction.'
    (args.folder/'movement_classification.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
