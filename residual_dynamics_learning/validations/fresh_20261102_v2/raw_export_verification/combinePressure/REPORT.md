# RL data report: combinePressure

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\RL_data_collection\20261102\combinePressure.mat`

## Coverage and quality

- Sensor rows: 6841, simulation interval 0.000..342.000 s.
- Valid orientation: 6828/6841; fresh paired frames: 6841/6841.
- Observable/valid theta: 6828/6841.
- Final command range: [0 0 0] to [2 2.027 2.027] bar.
- Sensor host elapsed: 347.662462 s; median interval 0.044909 s; maximum interval 0.375119 s.
- Overrun rows: 11.

## Experiment intervals

    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>trial_id</strong>    <strong>phase_id</strong>    <strong>channel_id</strong>    <strong>observation_rows</strong>
    <strong>__________________</strong>    <strong>_________________</strong>    <strong>________</strong>    <strong>________</strong>    <strong>__________</strong>    <strong>________________</strong>

               0                 4.95             0           0            0                100      
               5                 9.95             1           1            0                100      
              10                 72.8             1           5            0               1257      
           72.85                 77.8             1           3            0                100      
           77.85                 87.8             1           4            0                200      
           87.85                 92.8             2           1            0                100      
           92.85               155.65             2           5            0               1257      
           155.7               160.65             2           3            0                100      
           160.7               170.65             2           4            0                200      
           170.7               175.65             3           1            0                100      
           175.7               238.45             3           5            0               1256      
           238.5               243.45             3           3            0                100      
           243.5               253.45             3           4            0                200      
           253.5               258.45             4           1            0                100      
           258.5                321.3             4           5            0               1257      
          321.35                326.3             4           3            0                100      
          326.35                336.3             4           4            0                200      
          336.35                  342             0           9            0                114      



## Per-trial measured ranges

    <strong>trial_id</strong>    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>sensor_rows</strong>    <strong>valid_orientation_rows</strong>    <strong>phi_min_rad</strong>    <strong>phi_max_rad</strong>    <strong>channel_id</strong>    <strong>max_cycle_id</strong>    <strong>p1_min_bar</strong>    <strong>p1_max_bar</strong>    <strong>p2_min_bar</strong>    <strong>p2_max_bar</strong>    <strong>p3_min_bar</strong>    <strong>p3_max_bar</strong>
    <strong>________</strong>    <strong>__________________</strong>    <strong>_________________</strong>    <strong>___________</strong>    <strong>______________________</strong>    <strong>___________</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>____________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>

       1                  5                  87.8             1657                 1657              0.015163         0.1776          0              2              0                2          0                2          0                2  
       2              87.85                170.65             1657                 1655              0.018977        0.39383          0              2              0           1.8333          0           1.6667          0           1.6667  
       3              170.7                253.45             1656                 1650              0.019197        0.80204          0              2              0              1.5          0                2          0                2  
       4              253.5                 336.3             1657                 1652              0.017854        0.91138          0              2              0           1.8333          0            2.027          0            2.027  



Transitions crossing host acquisition gaps > 0.100 s are flagged, not silently used as uniform steps. The threshold is an export-quality setting, not a measured controller period. First/invalid transitions are also rejected.
There are 9 host acquisition gaps above that threshold; see acquisition_gaps.csv.

## Signal interpretation

Commands are pressure proxies, not chamber measurements. des_theta/des_phi are forward-model predictions. theta_reference_rad/phi_reference_rad, when present, are requested circle/helix angles, not observed targets. pressureDesired_bar is the inverse-model calculation and may drive the circle/helix feedforward path; final NI commands remain authoritative actions. Inferred lengths are not independent strain measurements. Corrected tip XYZ uses measured sensor positions plus orientation and assumed offsets. Native CSVs retain invalid/held source values; observations masks invalid theta/phi/inferred lengths. Use explicit per-feature validity; do not assume the whole row is valid for every target.

## Export and timing

observations.csv uses NDI simulation timestamps; commands.csv preserves the union of NI command log times. Native files preserve original signal times and values. Never treat the two clocks as interchangeable. Previous-sample alignment avoids interpolation through dropouts. See column_dictionary.csv for source BlockPath, signal index and constant-tail assumptions.

## Metadata and limitations

A recording is not a complete training dataset. Complete session/calibration/source settings/payload metadata in the supplied manifest. Split training/validation/test by whole trials or sessions. Do not estimate actual pressure lag solely from command-to-motion data. Source MAT was not modified.

## Warnings

- Constant-tail hold used for kill_switch
- Optional signal absent: direction
