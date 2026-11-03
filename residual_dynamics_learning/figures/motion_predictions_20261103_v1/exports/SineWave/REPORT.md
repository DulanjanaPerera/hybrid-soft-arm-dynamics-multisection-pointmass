# RL data report: SineWave

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\RL_data_collection\20261103\SineWave.mat`

## Coverage and quality

- Sensor rows: 5941, simulation interval 0.000..297.000 s.
- Valid orientation: 5931/5941; fresh paired frames: 5941/5941.
- Observable/valid theta: 5931/5941.
- Final command range: [0 0 0] to [3 3 3] bar.
- Sensor host elapsed: 300.822840 s; median interval 0.044607 s; maximum interval 0.210025 s.
- Overrun rows: 2.

## Experiment intervals

    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>trial_id</strong>    <strong>phase_id</strong>    <strong>channel_id</strong>    <strong>observation_rows</strong>
    <strong>__________________</strong>    <strong>_________________</strong>    <strong>________</strong>    <strong>________</strong>    <strong>__________</strong>    <strong>________________</strong>

              0                  4.95             0           0            0               100       
              5                 42.65             1           5            1               754       
           42.7                 52.65             1           4            1               200       
           52.7                 90.35             2           5            1               754       
           90.4                100.35             2           4            1               200       
          100.4                138.05             3           5            2               754       
          138.1                148.05             3           4            2               200       
          148.1                185.75             4           5            2               754       
          185.8                195.75             4           4            2               200       
          195.8                233.45             5           5            3               754       
          233.5                243.45             5           4            3               200       
          243.5                281.15             6           5            3               754       
          281.2                291.15             6           4            3               200       
          291.2                   297             0           9            0               117       



## Per-trial measured ranges

    <strong>trial_id</strong>    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>sensor_rows</strong>    <strong>valid_orientation_rows</strong>    <strong>phi_min_rad</strong>    <strong>phi_max_rad</strong>    <strong>channel_id</strong>    <strong>max_cycle_id</strong>    <strong>p1_min_bar</strong>    <strong>p1_max_bar</strong>    <strong>p2_min_bar</strong>    <strong>p2_max_bar</strong>    <strong>p3_min_bar</strong>    <strong>p3_max_bar</strong>
    <strong>________</strong>    <strong>__________________</strong>    <strong>_________________</strong>    <strong>___________</strong>    <strong>______________________</strong>    <strong>___________</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>____________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>

       1                  5                 52.65              954                 954               0.018733        1.8722           1              3              0             3             0             0             0             0     
       2               52.7                100.35              954                 954               0.082696         1.894           1              3              0             3             0             0             0             0     
       3              100.4                148.05              954                 954               0.070866         1.888           2              3              0             0             0             3             0             0     
       4              148.1                195.75              954                 954               0.077216        1.9041           2              3              0             0             0             3             0             0     
       5              195.8                243.45              954                 951               0.070612         2.449           3              3              0             0             0             0             0             3     
       6              243.5                291.15              954                 947               0.079793        2.4639           3              3              0             0             0             0             0             3     



Transitions crossing host acquisition gaps > 0.100 s are flagged, not silently used as uniform steps. The threshold is an export-quality setting, not a measured controller period. First/invalid transitions are also rejected.
There are 1 host acquisition gaps above that threshold; see acquisition_gaps.csv.

## Signal interpretation

Commands are pressure proxies, not chamber measurements. des_theta/des_phi are forward-model predictions. theta_reference_rad/phi_reference_rad, when present, are requested circle/helix angles, not observed targets. pressureDesired_bar is the inverse-model calculation and may drive the circle/helix feedforward path; final NI commands remain authoritative actions. Inferred lengths are not independent strain measurements. Corrected tip XYZ uses measured sensor positions plus orientation and assumed offsets. Native CSVs retain invalid/held source values; observations masks invalid theta/phi/inferred lengths. Use explicit per-feature validity; do not assume the whole row is valid for every target.

## Export and timing

observations.csv uses NDI simulation timestamps; commands.csv preserves the union of NI command log times. Native files preserve original signal times and values. Never treat the two clocks as interchangeable. Previous-sample alignment avoids interpolation through dropouts. See column_dictionary.csv for source BlockPath, signal index and constant-tail assumptions.

## Metadata and limitations

A recording is not a complete training dataset. Complete session/calibration/source settings/payload metadata in the supplied manifest. Split training/validation/test by whole trials or sessions. Do not estimate actual pressure lag solely from command-to-motion data. Source MAT was not modified.

## Warnings

- Constant-tail hold used for kill_switch
- Optional signal absent: direction
