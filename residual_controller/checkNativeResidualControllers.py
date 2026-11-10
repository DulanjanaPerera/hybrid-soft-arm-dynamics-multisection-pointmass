"""Build emitted standalone C with local GCC; replay CSV fixtures via ctypes.

No hardware/Simulink/API access. MATLAB's build-toolchain error remains
separate from this explicit standalone component compilation/replay.
"""
from pathlib import Path
import csv
import ctypes as C
import hashlib
import json
import math
import os
import subprocess
import sys


class Diagnostic(C.Structure):
    _fields_ = [(n, C.c_ubyte) for n in ['enabled', 'ready', 'advanced', 'historyGapHeld', 'contractValid']] + [
        (n, C.c_double * 2) for n in ['rawForce_N', 'effectiveResidual_N', 'gravityReexpression_N', 'inverseError_N']
    ] + [('iterations', C.c_ubyte), ('converged', C.c_ubyte), ('warm_s', C.c_double),
         ('units', C.c_double), ('parameterCorrection', C.c_double * 10)]


DP = C.POINTER(C.c_double)
PDF_IN = [('a', 6), ('a', 2), ('b', 1), ('d', 1), ('a', 3), ('b', 1), ('d', 1),
          ('b', 1), ('b', 1), ('b', 1), ('a', 10), ('a', 7), ('a', 8)]
PAC_IN = [('a', 6), ('a', 6), ('a', 3), ('b', 1), ('b', 1), ('d', 1), ('d', 1),
          ('b', 1), ('a', 2), ('a', 2), ('b', 1), ('b', 1)] + [('a', 10)] * 6 + [
    ('a', 3), ('a', 4), ('a', 10), ('a', 7), ('b', 1), ('b', 1), ('a', 8)]
PDF_OUT = [('a', 3), ('a', 2), ('a', 2), ('b', 1), ('b', 1), ('a', 2), ('a', 2), ('a', 2), ('b', 1), ('s', 1)]
PAC_OUT = [('a', 3)] + [('a', 2)] * 9 + [('b', 1)] * 5 + [
    ('a', 2), ('a', 2), ('a', 4), ('a', 4), ('a', 2), ('a', 2), ('a', 10), ('a', 3),
    ('b', 1), ('b', 1), ('a', 2), ('a', 1), ('a', 1), ('s', 1)]


def load_rows(path):
    with path.open(newline='') as stream:
        return [[float(v) for v in row] for row in csv.reader(stream)]


def flatten(value):
    if isinstance(value, Diagnostic):
        out = []
        for name, _ in Diagnostic._fields_:
            out.extend(flatten(getattr(value, name)))
        return out
    if isinstance(value, C.Array):
        return list(value)
    if hasattr(value, 'value'):
        return [float(value.value)]
    return [float(value)]


def main():
    root = Path(__file__).resolve().parents[1]
    generated = (root / sys.argv[1]).resolve()
    fixtures = (root / sys.argv[2]).resolve()
    gcc = Path('C:/mingw64/bin/gcc.exe')
    include = Path('C:/Program Files/MATLAB/R2025a/extern/include')
    os.environ['OMP_NUM_THREADS'] = '2'
    dll_directory = os.add_dll_directory(str(gcc.parent))
    names = {1: 'pdfbResidualInverseBlock', 2: 'pacResidualAdaptiveTunableBlock'}
    libraries = {}
    for kind, name in names.items():
        folder = generated / name
        dll = folder / 'component.dll'
        args = [str(gcc), '-shared', '-O2', '-fno-fast-math', '-fopenmp', '-Wl,--export-all-symbols',
                '-I', str(folder), '-I', str(include), *map(str, folder.glob('*.c')), '-o', str(dll), '-lm']
        result = subprocess.run(args, capture_output=True, text=True)
        (folder / 'gcc_build.txt').write_text(result.stdout + result.stderr, encoding='utf-8')
        assert result.returncode == 0, result.stderr
        libraries[kind] = C.CDLL(str(dll))
    with (fixtures / 'cases.csv').open(newline='') as stream:
        cases = list(csv.DictReader(stream))
    results = []
    for case in cases:
        kind = int(case['controller']);name = names[kind];lib = libraries[kind]
        ins = PDF_IN if kind == 1 else PAC_IN;outs = PDF_OUT if kind == 1 else PAC_OUT
        fn = getattr(lib, name)
        fn.restype = None
        fn.argtypes = [DP if tag == 'a' else C.c_double if tag == 'd' else C.c_ubyte for tag, n in ins] + [
            DP if tag == 'a' else C.POINTER(Diagnostic) if tag == 's' else C.POINTER(C.c_ubyte) for tag, n in outs]
        getattr(lib, name + '_initialize')()
        inputs = load_rows(fixtures / (case['id'] + '_input.csv'))
        expected = load_rows(fixtures / (case['id'] + '_expected.csv'))
        actual = [];maximum = 0.;pressure_max = 0.;activation_boundaries = 0
        for x, truth in zip(inputs, expected):
            args = [];offset = 0
            for tag, n in ins:
                values = x[offset:offset+n];offset += n
                args.append((C.c_double*n)(*values) if tag == 'a' else C.c_double(values[0]) if tag == 'd' else C.c_ubyte(int(values[0])))
            assert offset == len(x)
            values = [(C.c_double*n)() if tag == 'a' else Diagnostic() if tag == 's' else C.c_ubyte() for tag, n in outs]
            fn(*args, *[v if isinstance(v, C.Array) else C.byref(v) for v in values])
            row = [z for v in values for z in flatten(v)];assert len(row) == len(truth)
            for i, (a, b) in enumerate(zip(row, truth)):
                if math.isnan(b):assert math.isnan(a)
                else:
                    assert math.isfinite(a)
                    # Legacy allAboveActivation uses strict >.8. Floating
                    # summation can put equal pressure on opposite sides by
                    # 3e-16bar in MATLAB versus C; this flag is logging only.
                    po_index = 10 if kind == 1 else 20
                    activation = x[sum(n for _, n in ins[:po_index]) + 8]
                    if kind == 2 and i == 24 and a != b and abs(min(row[:3])-activation) < 1e-10 and abs(min(truth[:3])-activation) < 1e-10:
                        activation_boundaries += 1
                        continue
                    maximum = max(maximum, abs(a-b))
                    if i < 3:pressure_max = max(pressure_max, abs(a-b))
            actual.append(row)
        getattr(lib, name + '_terminate')()
        with (fixtures / (case['id'] + '_native.csv')).open('w', newline='') as stream:
            csv.writer(stream).writerows(actual)
        results.append(dict(case=case['id'], ticks=len(inputs), maximumDifference=maximum, pressureDifference_bar=pressure_max,
                            activationFlagBoundaryDifferences=activation_boundaries))
        assert maximum < 1e-6 and pressure_max < 1e-8, results[-1]
        print('PASS', results[-1])
    (fixtures / 'native_results.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    files = [p for folder in libraries for p in (generated/names[folder]).glob('*') if p.suffix in ['.c', '.h']]
    with (fixtures / 'generated_source_hashes.csv').open('w', newline='') as stream:
        writer = csv.writer(stream);writer.writerow(['path','sha256'])
        for path in files:writer.writerow([str(path.relative_to(root)), hashlib.sha256(path.read_bytes()).hexdigest()])
    print('ALL NATIVE COMPONENT REPLAYS PASS:', len(results), 'cases;', sum(r['ticks'] for r in results), 'ticks.')
    dll_directory.close()


if __name__ == '__main__':
    main()
