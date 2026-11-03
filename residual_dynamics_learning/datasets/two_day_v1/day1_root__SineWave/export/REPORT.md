# RL data report: SineWave

Source: `C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\RL_data_collection\SineWave.mat`

## Coverage and quality

- Sensor rows: 12721, simulation interval 0.000..636.000 s.
- Valid orientation: 12710/12721; fresh paired frames: 12721/12721.
- Observable/valid theta: 12710/12721.
- Final command range: [0 0 0] to [3 3 3] bar.
- Sensor host elapsed: 648.314467 s; median interval 0.056674 s; maximum interval 4.499000 s.
- Overrun rows: 17.

## Experiment intervals

    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>trial_id</strong>    <strong>phase_id</strong>    <strong>channel_id</strong>    <strong>observation_rows</strong>
    <strong>__________________</strong>    <strong>_________________</strong>    <strong>________</strong>    <strong>________</strong>    <strong>__________</strong>    <strong>________________</strong>

               0                 4.95             0           0            0                100      
               5                 99.2             1           5            1               1885      
           99.25                109.2             1           4            1                200      
          109.25               203.45             2           5            1               1885      
           203.5               213.45             2           4            1                200      
           213.5                307.7             3           5            2               1885      
          307.75                317.7             3           4            2                200      
          317.75               411.95             4           5            2               1885      
             412               421.95             4           4            2                200      
             422                516.2             5           5            3               1885      
          516.25                526.2             5           4            3                200      
          526.25               620.45             6           5            3               1885      
           620.5               630.45             6           4            3                200      
           630.5                  636             0           9            0                111      



## Per-trial measured ranges

    <strong>trial_id</strong>    <strong>first_sample_sim_s</strong>    <strong>last_sample_sim_s</strong>    <strong>sensor_rows</strong>    <strong>valid_orientation_rows</strong>    <strong>phi_min_rad</strong>    <strong>phi_max_rad</strong>    <strong>channel_id</strong>    <strong>max_cycle_id</strong>    <strong>p1_min_bar</strong>    <strong>p1_max_bar</strong>    <strong>p2_min_bar</strong>    <strong>p2_max_bar</strong>    <strong>p3_min_bar</strong>    <strong>p3_max_bar</strong>
    <strong>________</strong>    <strong>__________________</strong>    <strong>_________________</strong>    <strong>___________</strong>    <strong>______________________</strong>    <strong>___________</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>____________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>

       1                   5                109.2             2085                 2085              0.019435        1.9366           1              3              0             3             0             0             0             0     
       2              109.25               213.45             2085                 2085              0.061734        1.9502           1              3              0             3             0             0             0             0     
       3               213.5                317.7             2085                 2085              0.057435        1.9362           2              3              0             0             0             3             0             0     
       4              317.75               421.95             2085                 2085              0.091195        1.9591           2              3              0             0             0             3             0             0     
       5                 422                526.2             2085                 2079              0.065144        2.4902           3              3              0             0             0             0             0             3     
       6              526.25               630.45             2085                 2080              0.070345        2.5057           3              3              0             0             0             0             0             3     



Transitions crossing host acquisition gaps > 0.100 s are flagged, not silently used as uniform steps. The threshold is an export-quality setting, not a measured controller period. First/invalid transitions are also rejected.
There are 12 host acquisition gaps above that threshold; see acquisition_gaps.csv.

## Signal interpretation

Commands are pressure proxies, not chamber measurements. des_theta/des_phi are forward-model predictions. theta_reference_rad/phi_reference_rad, when present, are requested circle/helix angles, not observed targets. pressureDesired_bar is the inverse-model calculation and may drive the circle/helix feedforward path; final NI commands remain authoritative actions. Inferred lengths are not independent strain measurements. Corrected tip XYZ uses measured sensor positions plus orientation and assumed offsets. Native CSVs retain invalid/held source values; observations masks invalid theta/phi/inferred lengths. Use explicit per-feature validity; do not assume the whole row is valid for every target.

## Export and timing

observations.csv uses NDI simulation timestamps; commands.csv preserves the union of NI command log times. Native files preserve original signal times and values. Never treat the two clocks as interchangeable. Previous-sample alignment avoids interpolation through dropouts. See column_dictionary.csv for source BlockPath, signal index and constant-tail assumptions.

## Metadata and limitations

A recording is not a complete training dataset. Complete session/calibration/source settings/payload metadata in the supplied manifest. Split training/validation/test by whole trials or sessions. Do not estimate actual pressure lag solely from command-to-motion data. Source MAT was not modified.

## Warnings

- Constant-tail hold used for kill_switch
- Optional signal absent: direction
