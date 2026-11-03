# Capacity candidate model card

Architecture selected: gru128; representative seed113.
File: candidate_force_model.mat; SHA256 `8f2e0bb55ea9399b62abd5abe0b8f6dd15ab0470a400ad3e593e1edcd8f9f545`.
Outputs two generalized residual forces in N on RHS of the fixed nominal dynamics; not three pressures.
Use sgCell with cfg.widths and separate hidden states per layer. Existing normalizer unchanged.
Pressure-dependent residual requires bounded implicit inverse and immutable candidate-evaluation state.
Old TRAIN only; old/Nov3 selection; Nov2 evaluation only, all previously inspected.
GRU128 improves 4/10 fixed long examples relative to GRU64, while others worsen; these correlated examples are not independent statistical trials. Circle free RMSE worsens on Nov2 0.623->0.930mm and Nov3 1.005->1.366mm. It improves several sine/helix examples. The short-horizon winner is therefore a candidate for measured-state model augmentation, with mixed long-horizon accuracy. Retain the current GRU64 reference; do not infer improved physical tracking from short prediction scores.
This is an offline research candidate. Existing released GRU64 retained. No controller/hardware qualification.
