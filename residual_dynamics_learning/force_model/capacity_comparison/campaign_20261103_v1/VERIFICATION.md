# Verification

All nine internal scores reproduced; all nine external short-path tables checked.
All ten representative GRU64 one-step replays exactly reproduce previous fine-tune Q/R/anchors.
Original old sequences, segments and normalizer unchanged; original source/raw/prepared hashes pass.
Representative whole-record one-step and completed long paths independently recomputed: maximum discrepancy 4.88498130835e-15 mm.
Any incomplete long paths retained explicitly as failures, never scored as successful partial paths.
All Nov3 long windows are held-out trials. No interpolation.
All three inverse state-purity checks and208/208 prescribed-target bound checks pass; numerical force-equation errors reported separately per architecture.
Model hashes and original research release preserved. No hardware or controller qualification.

## Artifact review

All seven PNG figure bundles visually inspected: readable axes/legends, consistent architecture colors, no missing model traces. FIG/PDF copies retained.
An overview-only CSV text-column import error was repaired using explicit string imports and a render-only retry. No fits, predictions, numerical metrics or checkpoint choices changed; initial failure and successful retry remain in the evaluation diary.
