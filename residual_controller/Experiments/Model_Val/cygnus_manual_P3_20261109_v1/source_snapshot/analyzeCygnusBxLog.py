"""Offline Cygnus/Aurora BX:0001 hex log audit; no serial/device access."""
from pathlib import Path
import argparse
import collections
import csv
import hashlib
import json
import math
import re
import shutil
import struct


def crc16(data):
    value = 0
    for byte in data:
        value ^= byte
        for _ in range(8):
            value = (value >> 1) ^ 0xA001 if value & 1 else value >> 1
    return value


def decode(line):
    data = bytes.fromhex(line[2:])
    assert data[:2] == b'\xc4\xa5'
    size = struct.unpack_from('<H', data, 2)[0]
    assert len(data) == size + 8
    assert crc16(data[:4]) == struct.unpack_from('<H', data, 4)[0]
    body = data[6:-2]
    assert crc16(body) == struct.unpack_from('<H', data, len(data)-2)[0]
    offset = 1
    tools = {}
    for _ in range(body[0]):
        handle, status = body[offset:offset+2]
        offset += 2
        assert status in (1, 2, 4), status
        pose = None
        if status == 1:
            pose = struct.unpack_from('<8f', body, offset)
            offset += 32
            assert all(math.isfinite(v) for v in pose)
        port = frame = None
        if status != 4:
            port, frame = struct.unpack_from('<II', body, offset)
            offset += 8
        assert handle not in tools
        tools[handle] = dict(status=status, pose=pose, port=port, frame=frame)
    assert offset+2 == len(body)
    return tools, struct.unpack_from('<H', body, offset)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('log', type=Path)
    ap.add_argument('folder', type=Path)
    args = ap.parse_args()
    args.folder.mkdir(exist_ok=False, parents=True)
    blob = args.log.read_bytes()
    text = blob.decode('utf-8', errors='replace')
    frame_hz = float(re.search(r'<<Param.Tracking.Frame Frequency=([\d.]+)', text)[1])
    track_hz = float(re.search(r'<<Param.Tracking.Track Frequency=([\d.]+)', text)[1])
    session = 0
    packets = []
    last_command = ''
    for number, line in enumerate(text.splitlines(), 1):
        if line.startswith('>>'):
            last_command = line[2:]
            if last_command.startswith('TSTART'):
                session += 1
        if re.fullmatch(r'<<[0-9a-fA-F]+', line) and line[2:6].lower() == 'c4a5':
            assert last_command.startswith('BX:0001'), last_command
            tools, system = decode(line)
            assert set(tools) == {10, 11}
            packets.append(dict(session=session, line=number, tools=tools, system=system))
    main_session = max(p['session'] for p in packets)
    main_packets = [p for p in packets if p['session'] == main_session]
    frame0 = main_packets[0]['tools'][10]['frame']
    first_both = next(i for i, p in enumerate(main_packets)
                      if all(p['tools'][h]['status'] == 1 for h in (10, 11)))
    rows = []
    for i, packet in enumerate(main_packets):
        base, tip = [packet['tools'][h] for h in (10, 11)]
        assert base['frame'] == tip['frame']
        row = dict(index=i, log_line=packet['line'], elapsed_s=(base['frame']-frame0)/frame_hz,
                   frame=base['frame'], startup=int(i < first_both),
                   base_status=base['status'], tip_status=tip['status'],
                   base_port_hex=f"{base['port']:08X}", tip_port_hex=f"{tip['port']:08X}",
                   base_out_of_volume=int(bool(base['port'] & 0x40)),
                   tip_out_of_volume=int(bool(tip['port'] & 0x40)),
                   base_partial=int(bool(base['port'] & 0x80)),
                   tip_partial=int(bool(tip['port'] & 0x80)), system_status=packet['system'])
        for prefix, tool in [('base', base), ('tip', tip)]:
            values = tool['pose'] if tool['pose'] is not None else [float('nan')]*8
            for field, value in zip(['qw', 'qx', 'qy', 'qz', 'x_mm', 'y_mm', 'z_mm', 'error'], values):
                row[prefix+'_'+field] = value
        rows.append(row)
    runtime = rows[first_both:]
    missing = [row for row in runtime if row['tip_status'] == 2]
    episodes = []
    for row in missing:
        if not episodes or row['index'] != episodes[-1][-1]['index']+1:
            episodes.append([])
        episodes[-1].append(row)
    events = []
    for seq in episodes:
        next_row = rows[seq[-1]['index']+1]
        events.append(dict(start_s=seq[0]['elapsed_s'], last_missing_s=seq[-1]['elapsed_s'],
                           next_valid_s=next_row['elapsed_s'], missing_reports=len(seq),
                           reported_span_to_next_valid_s=next_row['elapsed_s']-seq[0]['elapsed_s'],
                           start_frame=seq[0]['frame'], last_missing_frame=seq[-1]['frame'],
                           first_log_line=seq[0]['log_line'], last_log_line=seq[-1]['log_line']))
    deltas = [b['frame']-a['frame'] for a, b in zip(rows, rows[1:])]
    summary = dict(source=str(args.log.resolve()), source_sha256=hashlib.sha256(blob).hexdigest(),
                   bx_replies=len(packets), crc_header_passed=len(packets), crc_payload_passed=len(packets),
                   main_tracking_session=main_session, main_tracking_replies=len(rows),
                   frame_counter_hz=frame_hz, transform_update_hz=track_hz, elapsed_s=rows[-1]['elapsed_s'],
                   time_origin='First BX reply in final tracking session; frame difference / reported frame frequency',
                   startup_reports=first_both, runtime_reports=len(runtime),
                   runtime_base_missing=sum(r['base_status'] == 2 for r in runtime),
                   runtime_tip_missing=len(missing), runtime_tip_missing_percent=100*len(missing)/len(runtime),
                   runtime_tip_missing_episodes=len(events), maximum_consecutive_missing=max(map(len, episodes), default=0),
                   missing_first_last_elapsed_s=([missing[0]['elapsed_s'], missing[-1]['elapsed_s']] if missing else []),
                   missing_base_valid_reports=sum(r['base_status'] == 1 for r in missing),
                   missing_out_of_volume_reports=sum(r['tip_out_of_volume'] for r in missing),
                   missing_partial_reports=sum(r['tip_partial'] for r in missing),
                   system_status_counts=dict(collections.Counter(p['system'] for p in main_packets)),
                   frame_increment_counts=dict(collections.Counter(deltas)),
                   sessions=[dict(session=s, replies=sum(p['session'] == s for p in packets))
                             for s in sorted(set(p['session'] for p in packets))],
                   pressure_markers_present=False)
    event_fields=['start_s','last_missing_s','next_valid_s','missing_reports',
                  'reported_span_to_next_valid_s','start_frame','last_missing_frame',
                  'first_log_line','last_log_line']
    for name, records in [('tracking.csv', rows), ('missing_tip_reports.csv', missing), ('missing_episodes.csv', events)]:
        with (args.folder/name).open('w', newline='') as f:
            fields=list(records[0]) if records else (list(rows[0]) if name=='missing_tip_reports.csv' else event_fields)
            writer = csv.DictWriter(f, fieldnames=fields);writer.writeheader();writer.writerows(records)
    (args.folder/'summary.json').write_text(json.dumps(summary, indent=2))
    shutil.copy2(args.log, args.folder/'source_log.txt')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
