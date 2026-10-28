# PAC tuning checkpoint: 2026-10-28

Checkpoint on codex/pac-autotuning-framework before residual dynamics research.
The working research path returns to the established manual-run workflow;
continuous online tuning is retained as experimental infrastructure, not a
validated unattended hardware system. No additional gain campaign is queued.

Included scope:

- Length-based target policy, fixed-raw smoothness objective, bounded proposals
  and manual campaign improvements; completed length/smoothness campaign ledgers,
  raw recordings, reproducible reviews, plots and source snapshots.
- Original measured-derivative filter restored and alternative filter retained
  for reproducibility; matched filter review did not support the alternative.
- Continuous-tuning runtime/worker, generated/manual model snapshots, session
  configurations, physical smoke data and failure reviews. Root model variants
  are saved experimental snapshots, not interchangeable approved configurations.
- Offline regression suites and documentation, including research rate note.

Latest physical R033 smoke completed three5s circles, no tracking limiting or
latched PAC faults, but rejected coverage88.44% and errors above targets. Eight
tip-validity dropouts in the200-sample score window need investigation. Worker
rejection-ID delay corrected offline; corrected behavior is not hardware-tested.
Earlier pacing, build-path, socket ownership and mixed-time validity problems
and their evidence are retained in VALIDATION_RESULTS and specialized reviews.
Runtime ownership locks/generated build caches are excluded. No hardware model
was started by the assistant. Existing manual control/PDFB/PAC behavior has not
been replaced by a residual model.

Unrelated RL collection edits (including a previously staged rename), paper
drafts, temporary manuscript scripts, exports and .original backups are outside
this checkpoint and remain in the workspace. No remote push is requested.

Next branch: codex/residual-dynamics-learning. Learn missing dynamics as an
augmentation of the model, validate predictions offline on held-out runs, then
compare original/augmented-model PDFB and PAC using manual experiments. This
is residual model learning, not a separate residual feedback controller.

Pre-commit validation: all nine offline suites passed in MATLAB R2025a:
validatePacTuneEvaluator, validatePacTuneManualCampaign, validatePacTuneOptimizer,
validatePacTuneOptimizerAdapter, validatePacSmoothMeasuredDerivatives,
validatePacTuneOnline, validatePacOnlineRuntime, validatePacOnlineValidity
(saved R02 prefix), and validatePacOnlineRejection (saved R033 score).
These checks do not establish physical robustness or unattended safety.
