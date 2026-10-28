function v=pacOnlineSettingsVector(scale,rate,lower,upper,initial,options,history, ...
    pressure,supervisor,lambda,Ks,enableAdaptation,enableCommon,freezeHysteresis)
%PACONLINESETTINGSVECTOR Audited fixed tuning order, exactly 81 doubles.
%#codegen
% Infinity stop time is encoded as -1 to permit strict finite packet checks.
stop=supervisor(:);if isinf(stop(7))&&stop(7)>0,stop(7)=-1;end
v=[scale(:);rate(:);lower(:);upper(:);initial(:);options(:);history(:); ...
    pressure(:);stop;lambda(:);Ks(:);double(enableAdaptation); ...
    double(enableCommon);double(freezeHysteresis)];
assert(numel(v)==81);
end
