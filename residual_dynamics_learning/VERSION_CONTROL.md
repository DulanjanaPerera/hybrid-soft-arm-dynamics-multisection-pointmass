# Residual research version-control checkpoint

Checkpoint prepared2026-11-03. Frozen model weights, numerical prediction
results, selection/evaluation tables, reports, source/fixtures/snapshots and
paper-candidate PNG/PDF/editable FIG figures are versioned. Full raw recordings,
large generated observation/command tables, prepared training datasets and
intermediate optimizer checkpoints remain LOCAL, preserved pending the user's
deferred migration outside OneDrive. No artifact was deleted or moved.

ARTIFACT_CHECKPOINT_MANIFEST.csv inventories residual/NDI artifacts and raw
recording inputs, with disposition, size and SHA256. It is a checkpoint
inventory, not a live automatically updated list; it excludes itself. Existing
campaign artifact_hashes/source manifests remain authoritative for each
experiment. The retained local files are still needed to rerun dataset-based
checks/retraining; a fresh Git checkout alone is not the full data archive.
Code/data paths in historical artifacts refer to the original workspace.

Git selection policy retains .m/.py/.md/.json/.txt, plots and compact CSV
results; excludes regenerated observations/commands/model_predictions,
row/window index tables and interpolation previews. MAT weights, protocols and
evaluation evidence are retained; checkpoint_*, *_sequences, *_prepared and
training/prepared data containers remain local. Generated code-build outputs
are omitted except explanatory reports/logs. Raw source MAT inputs are listed
but are not newly added by this checkpoint. Previously tracked input files
remain tracked. Papers, tmp/paper_draft and operator SLX changes are preserved
outside these commits, including the existing staged SLX rename.

The nested .gitattributes preserves exact archived bytes so SHA256 provenance
survives checkout without newline conversion in residual/NDI namespaces. Other
historical canonical sources retain repository newline conventions; compare
archived source snapshots when byte hashes differ solely by line endings.

Current model paths/hashes and next controller sequence are in
project_docs/RESIDUAL_CONTROLLER_HANDOFF.md. Primary candidate GRU128 improves
short prediction, while GRU64 remains the reference for better circle free
runs. HIL update automation is deferred until after the initial controller
implementation/comparison. No hardware or Simulink execution or remote push.
