function [theta, phi] = commonPressureFeedbackReferenceBlock(t)
%#codegen
% Paste into the existing thetaCircle MATLAB Function block in a COPY of
% the feedback-linearization model. Connect theta and phi to the desired
% angle Mux and NDI comparison inputs; bypass the old circle/ramp sources.
% Stop Time: 240 s. Both reference steps are only 0.01 rad.

theta = 0.0;
phi = 0.0;
if ~isfinite(t) || t < 0.0
    return;
end

if t >= 90.0 && t < 130.0
    theta = 0.0; % U-direction bend.
    if t < 95.0
        u = (t-90.0)/5.0;
        phi = 0.01*u*u*(3.0-2.0*u);
    elseif t < 125.0
        phi = 0.01;
    else
        u = (t-125.0)/5.0;
        phi = 0.01*(1.0-u*u*(3.0-2.0*u));
    end
elseif t >= 160.0 && t < 200.0
    theta = pi/2; % V-direction bend.
    if t < 165.0
        u = (t-160.0)/5.0;
        phi = 0.01*u*u*(3.0-2.0*u);
    elseif t < 195.0
        phi = 0.01;
    else
        u = (t-195.0)/5.0;
        phi = 0.01*(1.0-u*u*(3.0-2.0*u));
    end
end
end
