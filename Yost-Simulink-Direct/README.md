# Yost direct Simulink reader

Windows desktop Simulink (Normal mode), one dedicated wireless dongle per sensor.
No UDP bridge or Java runtime is required. Source is independent of the GUI project.

In MATLAB, change to this folder and run `setup_yost`. This generates the library
and demo .slx files, runs synthetic checks, and opens the demo. A precompiled
R2025a Windows MEX is included. To rebuild with your selected MATLAB C++ compiler,
run `build_yost` (`mex -setup C++` first if necessary).
The demo starts in synthetic-data mode; double-click the
Yost block to edit its configuration and turn synthetic mode off for hardware.
Close the UDP bridge and Yost Sensor Suite before opening their COM ports here.
COM28 is the dongle; COM27 is currently used only for charging.

The mask accepts a numeric sensor table, one row per sensor:
`[COM_number Yost_ID output_ID requested_Hz channel_mask magnetometer]`.
Example: `[28 1 1 100 3 0; 29 1 2 100 3 0]` (replace COM29 with your second dongle).
Yost IDs are 0..14, output IDs 0..255, rate is 1..1000 Hz or 0 for maximum.
Channel mask bits: quaternion=1, angular velocity=2, corrected acceleration=4,
linear acceleration=8. Sum the selected bits; 3 means quaternion plus gyro.
Magnetometer is 0/1. Gyroscope and accelerometer fusion are enabled; current
orientation filter and existing tare/mount offsets are preserved. No flash commit.
Port enumeration: `yost_list_ports` (does not open ports).

Each output matrix has one COLUMN per configured sensor, in table order:

| Port | Size | Contents |
|---|---|---|
| Quaternion | 4 x N | x,y,z,w (scalar last), existing sensor tare applied |
| Gyro | 3 x N | x,y,z in rad/s, command 0x21 |
| Accel | 3 x N | corrected x,y,z in g, command 0x27 |
| LinearAccel | 3 x N | gravity removed, global frame, g, command 0x29 |
| Info | 9 x N | output ID; sensor timestamp (microseconds); new sample flag; fresh flag; host sample age (seconds); accepted sample count; samples superseded before block read; parser discarded bytes; reader state |

Reader state: 1=waiting, 2=receiving, -1=serial/worker error. Fresh requires a
received sample, healthy reader, and age below the mask stale-time limit.
All ports use double. Before the first sample and for disabled channels, data
are NaN. Last values are held during a gap: consumers MUST use fresh/new flags.
Attach Terminator blocks to unwanted ports; the channel mask reduces radio traffic.
The device timestamp is raw uint32 represented exactly as double, wraps about
every 71.6 minutes and is not synchronized between dongles. Use host age for freshness.

Model sample time and requested IMU rate are separate. Default block time is
0.01 s and device rate is 100 Hz. The reader stores the latest sample, not every
sample for logging. Superseded samples are not proof of wireless loss. Use separate
logging if every acquired frame is required. Desktop simulation time is not wall
time: use simulation pacing (enabled at 1x in the example); Windows is not a hard
real-time scheduler. Synthetic mode advances on simulation time and opens no COM.

On simulation start the block configures each dongle, then starts a serial worker
per dongle. Model outputs only copy a snapshot. Stop joins workers and releases
COM handles. Start failures close already opened ports. Initialization can take
seconds while configuring hardware. Do not run the same COM in multiple blocks.
Fast Restart, Accelerator, Rapid Accelerator, standalone code generation, Simulink
Desktop Real-Time and Simulink Real-Time are not supported by this initial package.

The serial parser checks header shape, logical ID, exact selected payload size and
finite floats, and attempts byte resynchronization. It does not provide a CRC or
guaranteed recovery from arbitrary corruption. NDI proximity/radio gaps can persist.
This is a prototype pending live hardware validation; a passing synthetic test is
not evidence of correct wireless behavior or timing. Calibration/tare controls
remain available in the separate GUI; stop the model before using that GUI.

Files: `yost_reader.hpp` (Windows transport/parser), `yost_imu_sfun.cpp` (Simulink
adapter), `build_yost.m`, `create_yost_models.m`, `validate_yost.m` (synthetic checks),
`setup_yost.m` (entry point), `yost_list_ports.m`, `build_native.ps1`, and
`test_parser.cpp` (standalone parser checks).

## Validation status, 2026-11-09

- C++ syntax check against installed R2025a headers: passed.
- Standalone MEX compilation/linking with MinGW 13.1 and R2025a libraries: passed.
- Native parser tests: passed 521 channel-mask/frame-split cases, noise/wrong-ID/
  NaN rejection and recovery, timestamp rollover, and 10,000 batched frames.
- MATLAB batch startup failed with `System Error: File system inconsistency`
  before executing any MATLAB code. Therefore MEX loading, generated models,
  and `validate_yost` have NOT yet been verified inside Simulink. The .slx files
  will be generated when you run `setup_yost` in your working MATLAB session.
- No physical COM port was opened; live sensor startup, release and timing are
  still to be tested. The existing UDP project and MATLAB sessions were untouched.

Reference: [MathWorks simulation pacing options](https://www.mathworks.com/help/simulink/slref/simulationpacingoptions.html).
