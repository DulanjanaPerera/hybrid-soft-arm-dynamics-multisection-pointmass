# RL data report: combinePressure

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\RL_data_collection\20261103\combinePressure.mat`

## Coverage and quality

- Sensor rows: 3821, simulation interval 0.000..191.000 s.
- Valid orientation: 3811/3821; fresh paired frames: 3821/3821.
- Observable/valid theta: 3811/3821.
- Final command range: [0 0 0] to [3 3 3] bar.
- Sensor host elapsed: 193.808955 s; median interval 0.044746 s; maximum interval 0.256593 s.
- Overrun rows: 2.

## Experiment intervals

    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>trial_id</strong>    <strong>phase_id</strong>    <strong>channel_id</strong>    <strong>observation_rows</strong>
    <strong>__________________</strong>    <strong>_________________</strong>    <strong>________</strong>    <strong>________</strong>    <strong>__________</strong>    <strong>________________</strong>

               0                 4.95             0           0            0               100       
               5                 9.95             1           1            0               100       
              10                 35.1             1           5            0               503       
           35.15                 40.1             1           3            0               100       
           40.15                 50.1             1           4            0               200       
           50.15                 55.1             2           1            0               100       
           55.15                80.25             2           5            0               503       
            80.3                85.25             2           3            0               100       
            85.3                95.25             2           4            0               200       
            95.3               100.25             3           1            0               100       
           100.3               125.35             3           5            0               502       
           125.4               130.35             3           3            0               100       
           130.4               140.35             3           4            0               200       
           140.4               145.35             4           1            0               100       
           145.4                170.5             4           5            0               503       
          170.55                175.5             4           3            0               100       
          175.55                185.5             4           4            0               200       
          185.55                  191             0           9            0               110       



## Per-trial measured ranges

    <strong>trial_id</strong>    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>sensor_rows</strong>    <strong>valid_orientation_rows</strong>    <strong>phi_min_rad</strong>    <strong>phi_max_rad</strong>    <strong>channel_id</strong>    <strong>max_cycle_id</strong>    <strong>p1_min_bar</strong>    <strong>p1_max_bar</strong>    <strong>p2_min_bar</strong>    <strong>p2_max_bar</strong>    <strong>p3_min_bar</strong>    <strong>p3_max_bar</strong>
    <strong>________</strong>    <strong>__________________</strong>    <strong>_________________</strong>    <strong>___________</strong>    <strong>______________________</strong>    <strong>___________</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>____________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>

       1                  5                  50.1              903                 903               0.024036        0.25317          0              2              0                3          0                3          0                3  
       2              50.15                 95.25              903                 903               0.024416        0.42784          0              2              0           1.8333          0           1.6667          0           1.6667  
       3               95.3                140.35              902                 897               0.029456        0.85984          0              2              0              1.5          0                2          0                2  
       4              140.4                 185.5              903                 898               0.027421        0.96777          0              2              0           1.8333          0            2.027          0            2.027  



Transitions crossing host acquisition gaps > 0.100 s are flagged, not silently used as uniform steps. The threshold is an export-quality setting, not a measured controller period. First/invalid transitions are also rejected.
There are 2 host acquisition gaps above that threshold; see acquisition_gaps.csv.

## Signal interpretation

Commands are pressure proxies, not chamber measurements. des_theta/des_phi are forward-model predictions. theta_reference_rad/phi_reference_rad, when present, are requested circle/helix angles, not observed targets. pressureDesired_bar is the inverse-model calculation and may drive the circle/helix feedforward path; final NI commands remain authoritative actions. Inferred lengths are not independent strain measurements. Corrected tip XYZ uses measured sensor positions plus orientation and assumed offsets. Native CSVs retain invalid/held source values; observations masks invalid theta/phi/inferred lengths. Use explicit per-feature validity; do not assume the whole row is valid for every target.

## Export and timing

observations.csv uses NDI simulation timestamps; commands.csv preserves the union of NI command log times. Native files preserve original signal times and values. Never treat the two clocks as interchangeable. Previous-sample alignment avoids interpolation through dropouts. See column_dictionary.csv for source BlockPath, signal index and constant-tail assumptions.

## Metadata and limitations

A recording is not a complete training dataset. Complete session/calibration/source settings/payload metadata in the supplied manifest. Split training/validation/test by whole trials or sessions. Do not estimate actual pressure lag solely from command-to-motion data. Source MAT was not modified.

## Warnings

- Constant-tail hold used for kill_switch
- Optional signal absent: direction
