# November3 exports and NDI gap audit

All five raw recordings remain in RL_data_collection/20261103 unchanged.
Exports under this folder are for quality audit, not an assigned training split.
November2 motion plots remain in the adjacent motion_predictions_20261102_v1
folder. Additional November3 frozen-model plots are now complete under
frozen_model_predictions_v1/REPORT.md. No November3 fitting occurred; data have
now been model-evaluated/visualized and still have no assigned training split.

Authoritative NDI episode/sample counts are in short_gap_audit/gap_candidates.csv
and REPORT.md:55episodes/57samples.53 single-frame previews qualify for<=.12s
valid-to-valid span. Neither original flags nor model pipelines were changed.
The initial preparation console accidentally printed a reused cycle-count
variable as gap count; that print statement is corrected in source. Exported
observations and gap episode tables were unaffected. The old ndi_gaps.csv
optional reason column is unavailable; actual tip status and reject MASK values
are preserved in short_gap_audit/gap_candidates.csv and observations.csv.

Recorded settings are in recorded_trials.csv, not inferred from current source
block defaults: circles10s, requested phi[.2,.8,1,1.5,2.5]; sine cycle spacing
12.575s; helix10s/turn, requested peak phi2.5. These are not2s-circle qualification.

Interpolation preview/proxy results are not measured truth at real dropouts and
not evidence of learning improvement. Declare train/heldout roles before any
later fitting; retain raw source/calibration/day metadata.
