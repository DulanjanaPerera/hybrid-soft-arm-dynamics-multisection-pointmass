function [position,valid] = ndiValidSample(position,valid,frame,previousFrame)
% Reject missing, nonfinite, repeated, or mismatched frames before geometry.
valid = logical(valid) & all(isfinite(position),1) & isfinite(frame) ...
    & (frame ~= previousFrame);
if frame(1) ~= frame(2)
    valid(:) = false;
end
position(:,~valid) = NaN;
end
