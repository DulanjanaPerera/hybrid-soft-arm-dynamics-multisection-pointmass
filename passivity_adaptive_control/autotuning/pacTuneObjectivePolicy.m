function [metric,normalizer]=pacTuneObjectivePolicy(cfg)
%PACTUNEOBJECTIVEPOLICY Objective is immutable per campaign; old configs use q.
metric='qVectorRmse_mm';
if isfield(cfg,'objectiveMetric'),metric=char(cfg.objectiveMetric);end
switch metric
    case 'qVectorRmse_mm'
        normalizer=cfg.lengthNormalizer_mm;
    case 'velocityErrorRateRmse_mPerS2'
        assert(isfield(cfg,'smoothnessNormalizer_mPerS2'),'pacTune:Objective','Missing smoothness normalizer.');
        normalizer=cfg.smoothnessNormalizer_mPerS2;
    otherwise
        error('pacTune:Objective','Unsupported objective metric: %s.',metric);
end
assert(isnumeric(normalizer)&&isscalar(normalizer)&&isfinite(normalizer)&&normalizer>0, ...
    'pacTune:Objective','Objective normalizer must be positive and finite.');
end
