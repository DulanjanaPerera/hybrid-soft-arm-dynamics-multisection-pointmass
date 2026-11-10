"""Independently compile/replay emitted passive-block C. No model/device I/O."""
from pathlib import Path
import csv
import ctypes as C
import json
import math
import os
import subprocess
import sys


class Context(C.Structure):
    _fields_ = [('valid', C.c_ubyte), ('measured', C.c_double * 4),
                ('historyPre', C.c_double * 128), ('pressure_bar', C.c_double * 3)] + [
        (n, C.c_double) for n in ['featureDt_s', 'horizon_s', 'sampleTime_host_s', 'units']]


class Diagnostic(C.Structure):
    _fields_ = [(n, C.c_ubyte) for n in ['valid', 'networkEvaluations',
        'allocationLimited', 'slewLimited', 'nominalAllocationLimited',
        'nominalSlewLimited', 'frozenBalanceSatisfied']] + [
        (n, C.c_double * 2) for n in ['targetForce_N', 'residualForce_N',
        'correctedForceRequest_N', 'frozenTotalForce_N', 'allocationError_N',
        'slewError_N', 'frozenBalanceError_N', 'nominalBalanceError_N']] + [
        (n, C.c_double * 3) for n in ['residualInputPressure_bar',
        'correctedPressureBeforeSlew_bar', 'nominalPressureBeforeSlew_bar',
        'pressureDifferenceBeforeSlew_bar', 'pressureDifference_bar']] + [
        ('sampleTime_host_s', C.c_double), ('horizon_s', C.c_double)]


def rows(path):
    with path.open(newline='') as f:
        return [[float(x) for x in row] for row in csv.reader(f)]


def flat(s):
    return [float(x) for name, kind in s._fields_
            for x in (list(getattr(s, name)) if issubclass(kind, C.Array)
                      else [getattr(s, name)])]


def main():
    root = Path(__file__).resolve().parents[1]
    generated, fixtures = [(root / p).resolve() for p in sys.argv[1:3]]
    name = sys.argv[3] if len(sys.argv) > 3 else 'rlResidualSinglePassShadowBlock'
    assert name in ['rlResidualSinglePassShadowBlock', 'rlResidualSinglePassCollectionShadowBlock',
                    'rlResidualDesiredPressureBlock']
    reference_size = 6 if name == 'rlResidualDesiredPressureBlock' else 2
    gcc = Path('C:/mingw64/bin/gcc.exe')
    dll_path = os.add_dll_directory(str(gcc.parent))
    command = [str(gcc), '-shared', '-O2', '-fno-fast-math', '-Wl,--export-all-symbols',
               '-I', str(generated), '-I', 'C:/Program Files/MATLAB/R2025a/extern/include',
               *[str(p) for p in generated.glob('*.c')],
               '-o', str(generated / 'single_pass.dll'), '-lm']
    build = subprocess.run(command, capture_output=True, text=True)
    (fixtures / 'gcc_build.txt').write_text(build.stdout + build.stderr, encoding='utf-8')
    assert build.returncode == 0, build.stderr
    library = C.CDLL(str(generated / 'single_pass.dll'))
    entry = getattr(library, name)
    dp = C.POINTER(C.c_double)
    entry.argtypes = [C.POINTER(Context), dp, dp, C.c_ubyte, dp, dp, dp, C.POINTER(Diagnostic)]
    entry.restype = None
    getattr(library, name + '_initialize')()
    inputs = rows(fixtures / 'inputs.csv'); expected = rows(fixtures / 'expected.csv')
    assert len(inputs) == len(expected) and inputs
    maximum = 0.; actual = []
    for i, (x, reference) in enumerate(zip(inputs, expected)):
        ctx = Context(); offset = 0
        for name, kind in Context._fields_:
            if issubclass(kind, C.Array):
                size = kind._length_; setattr(ctx, name, kind(*x[offset:offset+size])); offset += size
            else:
                setattr(ctx, name, int(x[offset]) if kind == C.c_ubyte else x[offset]); offset += 1
        p, pn, diagnostic = (C.c_double*3)(), (C.c_double*3)(), Diagnostic()
        options_offset = offset + reference_size
        assert len(x) == options_offset + 13
        entry(C.byref(ctx), (C.c_double*reference_size)(*x[offset:options_offset]),
              (C.c_double*10)(*x[options_offset:options_offset+10]), int(x[options_offset+10]),
              (C.c_double*2)(*x[options_offset+11:options_offset+13]), p, pn, C.byref(diagnostic))
        result = list(p) + list(pn) + flat(diagnostic)
        assert len(result) == len(reference)
        for j, (got, want) in enumerate(zip(result, reference)):
            if math.isnan(want):
                assert math.isnan(got), (i, j, got, want)
            else:
                difference = abs(got-want)
                assert difference <= 1e-8*max(1., abs(want)), (i, j, got, want)
                maximum = max(maximum, difference)
        actual.append(result)
    with (fixtures / 'native_actual.csv').open('w', newline='') as f:
        csv.writer(f).writerows(actual)
    result = dict(passed=True, rows=len(actual), maximumAbsoluteDifference=maximum,
                  scope='Offline native component replay, no SLX/hardware/full-loop qualification')
    (fixtures / 'native_result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result))
    del dll_path


if __name__ == '__main__':
    main()
