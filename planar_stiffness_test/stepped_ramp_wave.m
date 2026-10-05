function [p1,p2,p3] = stepped_ramp_wave(t)
%#codegen
% Three P1 load/unload cycles; commands in bar. Total duration 400 s.
levels = [0 .5 1 1.5 2 2.5 3 2.5 2 1.5 1 .5 0];
hold_s = 10;
cycles = 3;
p1=0; p2=0; p3=0;
if ~isfinite(t) || t<0 || t>=hold_s*numel(levels)*cycles, return; end
slot=floor(t/hold_s);
p1=levels(mod(slot,numel(levels))+1);
end
