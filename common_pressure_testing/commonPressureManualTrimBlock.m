function [P1_bar, P2_bar, P3_bar, phaseID, limitFlag] = ...
    commonPressureManualTrimBlock(s_bar, a_bar, trialID, enablePressure)
%#codegen
% Paste into a MATLAB Function block.
%
% Inputs are tunable signals, not edits to this source during a run:
%   s_bar, a_bar: absolute differential settings in bar.
%   trialID: manually chosen label for the current trial.
%   enablePressure: 0 outputs zero; 1 outputs active pressure targets.
%
% The active mean command is fixed at 1.6 bar:
%   P1 = c-2*s/3, P2 = c+s/3-a, P3 = c+s/3+a.
% s=a=-0.15 gives the observed near-target [1.7,1.7,1.4] bar.
% Apply a Simulink Rate Limiter AFTER the three pressures are muxed;
% log its output as the delivered software command. This block itself
% outputs targets and has no time-dependent ramp or NDI feedback.

common_bar = 1.6;
activeFloor_bar = 1.2; % Provisional tested active-region floor.
ceiling_bar = 3.0;

P1_bar = 0.0;
P2_bar = 0.0;
P3_bar = 0.0;
phaseID = 0.0;
limitFlag = 0.0;

if ~isfinite(s_bar) || ~isfinite(a_bar) || ...
        ~isfinite(trialID) || ~isfinite(enablePressure)
    phaseID = -1.0;
    limitFlag = 1.0;
    return;
end

phaseID = round(trialID);
if enablePressure <= 0.5
    return;
end

wanted = [common_bar - (2.0/3.0)*s_bar, ...
          common_bar + (1.0/3.0)*s_bar - a_bar, ...
          common_bar + (1.0/3.0)*s_bar + a_bar];
limited = min(max(wanted, activeFloor_bar), ceiling_bar);
if any(abs(limited - wanted) > 1.0e-12)
    limitFlag = 1.0; % Discard a clipped trial from trim estimation.
end

P1_bar = limited(1);
P2_bar = limited(2);
P3_bar = limited(3);
end
