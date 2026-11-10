"""Build and verify emitted shadow-block C offline. No Simulink/device access."""
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


class ForwardDiagnostic(C.Structure):
    _fields_ = [(n, C.c_ubyte) for n in ['status', 'accepted', 'advanced', 'ready',
                'historyReset', 'predictionValid', 'scoreValid', 'retimedScoreValid', 'inverseContextValid']] + [
        (n, C.c_double * 2) for n in ['qMeasured_m', 'velocity_mPerS', 'rawForce_N']] + [
        (n, C.c_double * 6) for n in ['errorAtMeasurement_m', 'retimedError_m']] + [
        (n, C.c_double) for n in ['sampleTime_host_s', 'predictionTarget_host_s', 'scoreOrigin_host_s',
                                'dtPhysical_s', 'commandAge_sim_s', 'horizonMismatch_s']] + [
        ('pressureChangeSinceOrigin_bar', C.c_double * 3)] + [
        (n, C.c_double) for n in ['velocityWarm_s', 'neuralWarm_s', 'maximumTrainingZ', 'units']] + [
        ('segment', C.c_uint)]


class InverseDiagnostic(C.Structure):
    _fields_ = [(n, C.c_ubyte) for n in ['valid', 'converged', 'limited', 'nominalLimited', 'iterations']] + [
        (n, C.c_double * 2) for n in ['targetForce_N', 'predictedForce_N', 'inverseError_N', 'nominalInverseError_N']] + [
        ('pressureDifference_bar', C.c_double * 3), ('sampleTime_host_s', C.c_double)]


def flat(value):
    if isinstance(value, C.Structure):
        return [v for name, _ in value._fields_ for v in flat(getattr(value, name))]
    if isinstance(value, C.Array):
        return list(value)
    return [float(value)]


def rows(path):
    with path.open(newline='') as f:
        return [[float(v) for v in line] for line in csv.reader(f)]


def main():
    root = Path(__file__).resolve().parents[1]
    generated, fixture = [(root / p).resolve() for p in sys.argv[1:3]]
    gcc = Path('C:/mingw64/bin/gcc.exe')
    dll_directory = os.add_dll_directory(str(gcc.parent))  # retain handle until replay ends
    libraries = {}
    for kind, entry in [('forward', 'rlResidualShadowBlock'), ('inverse', 'rlResidualInverseShadowBlock')]:
        folder = generated / entry
        command = [str(gcc), '-shared', '-O2', '-fno-fast-math', '-Wl,--export-all-symbols',
                   '-I', str(folder), '-I', 'C:/Program Files/MATLAB/R2025a/extern/include',
                   *[str(p) for p in folder.glob('*.c')], '-o', str(folder / 'component.dll'), '-lm']
        build = subprocess.run(command, capture_output=True, text=True)
        (folder / 'gcc_build.txt').write_text(build.stdout + build.stderr, encoding='utf-8')
        assert build.returncode == 0, build.stderr
        libraries[kind] = C.CDLL(str(folder / 'component.dll'))
    dp = C.POINTER(C.c_double)
    forward = libraries['forward'].rlResidualShadowBlock
    forward.restype = None
    forward.argtypes = [dp, C.c_ubyte, C.c_ubyte, C.c_double, C.c_double, dp,
                       C.c_double, C.c_double, C.c_ubyte, dp, dp, dp, dp,
                       C.POINTER(ForwardDiagnostic), C.POINTER(Context)]
    inverse = libraries['inverse'].rlResidualInverseShadowBlock
    inverse.restype = None
    inverse.argtypes = [C.POINTER(Context), dp, dp, C.c_ubyte, dp, dp, dp, C.POINTER(InverseDiagnostic)]
    results = []
    with (fixture / 'cases.csv').open(newline='') as f:
        cases = list(csv.DictReader(f))
    for case in cases:
        inputs = rows(fixture / (case['id'] + '_inputs.csv'))
        expected = rows(fixture / (case['id'] + '_expected.csv'))
        if case['kind'] == 'forward':
            libraries['forward'].rlResidualShadowBlock_initialize()
        else:
            libraries['inverse'].rlResidualInverseShadowBlock_initialize()
        actual = []
        maximum = 0.
        for x, reference in zip(inputs, expected):
            if case['kind'] == 'forward':
                a, b, c = [(C.c_double * 6)() for _ in range(3)]
                d, context = ForwardDiagnostic(), Context()
                forward((C.c_double*3)(*x[:3]), int(x[3]), int(x[4]), x[5], x[6],
                        (C.c_double*3)(*x[7:10]), x[10], x[11], int(x[12]),
                        (C.c_double*11)(*x[13:24]), a, b, c, C.byref(d), C.byref(context))
                result = list(a) + list(b) + list(c) + flat(d) + flat(context)
            else:
                context = Context()
                offset = 0
                for name, typ in Context._fields_:
                    if issubclass(typ, C.Array):
                        n = typ._length_
                        setattr(context, name, typ(*x[offset:offset+n])); offset += n
                    else:
                        setattr(context, name, int(x[offset]) if typ == C.c_ubyte else x[offset]); offset += 1
                a, b, d = (C.c_double*3)(), (C.c_double*3)(), InverseDiagnostic()
                inverse(C.byref(context), (C.c_double*2)(*x[offset:offset+2]),
                        (C.c_double*10)(*x[offset+2:offset+12]), int(x[offset+12]),
                        (C.c_double*8)(*x[offset+13:offset+21]), a, b, C.byref(d))
                result = list(a) + list(b) + flat(d)
            assert len(result) == len(reference)
            for j, (got, want) in enumerate(zip(result, reference)):
                if math.isnan(want):
                    assert math.isnan(got), (case['id'], len(actual), j, got, want)
                else:
                    difference = abs(got-want)
                    assert difference <= 1e-8*max(1., abs(want)), (case['id'], len(actual), j, got, want)
                    maximum = max(maximum, difference)
            actual.append(result)
        with (fixture / (case['id'] + '_native.csv')).open('w', newline='') as f:
            csv.writer(f).writerows(actual)
        results.append(dict(case=case['id'], rows=len(actual), maxAbsoluteDifference=maximum, passed=True))
    (fixture / 'native_results.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    print(json.dumps(results, indent=2))
    dll_directory.close()


if __name__ == '__main__':
    main()
