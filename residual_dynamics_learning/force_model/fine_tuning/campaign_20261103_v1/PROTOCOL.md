# Predeclared fine-tune experiment

Exact GRU64 weight warm start: 8150352dd47230a2a4094b240e1e886b641b8c22ebad0cbfe300f63e9f526c6a
Old October19/27 splits and normalizer unchanged. Nov3 middle whole trial per channel reserved, no interpolation. Nov2 never used for fitting/selection.
Adam restarted, learning rate 0.0001,1200 extra updates per arm; seed113. Mixed16 old/16 new streams; old-only32 old. Same loss/history/physics.
Selection=.5old+.5Nov3, each .7oneStep+.3halfSecond macro RMSE. Old guard<=1.10 initial old score. Initial weights remain eligible.
Initial old/new/composite: 0.544112625/0.87840448/0.711258553 mm.
All models frozen before Nov2 reading. No hardware/controller change.
