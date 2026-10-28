function [ready,required,active]=pacTuneTargetPolicy(cfg)
%PACTUNETARGETPOLICY Versioned acceptance policy; NaN disables optional limits.
if strcmp(cfg.schemaVersion,'pac-tune-1')
    required={'qVectorRmse_mm','thetaMae_rad','phiMae_rad','velocityRmse_mPerS','limitedFraction'};
    allowed=required;
elseif strcmp(cfg.schemaVersion,'pac-tune-2')
    required={'q2Rmse_mm','q3Rmse_mm','velocityRmse_mPerS', ...
        'velocityErrorRateRmse_mPerS2','limitedFraction'};
    allowed=[required {'qVectorRmse_mm','thetaMae_rad','phiMae_rad','thetaRateStd_radPerS'}];
else
    error('pacTune:Schema','Unknown evaluator target schema.');
end
assert(all(isfield(cfg.targets,required)),'pacTune:Targets','Missing required target fields.');
active={};
for f=fieldnames(cfg.targets)'
    x=cfg.targets.(f{1});
    assert(ismember(f{1},allowed)&&isnumeric(x)&&isscalar(x)&& ...
        (isnan(x)||(isfinite(x)&&x>=0)),'pacTune:Targets', ...
        'Target %s must be an allowed scalar nonnegative limit or NaN.',f{1});
    if isfinite(x),active{end+1}=f{1};end %#ok<AGROW>
end
if strcmp(cfg.schemaVersion,'pac-tune-2')&&isfinite(cfg.targets.limitedFraction)
    assert(cfg.targets.limitedFraction<=1,'pacTune:Targets','limitedFraction must lie in [0,1].');
end
ready=all(cellfun(@(f)isfinite(cfg.targets.(f)),required));
end
