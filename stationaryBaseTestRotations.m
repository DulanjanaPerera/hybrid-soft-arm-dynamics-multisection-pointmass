function rotations = stationaryBaseTestRotations()
% Identity, gravity-axis yaw, signed pitch/roll, and fixed sampled rotations.
rx=@(a)[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
ry=@(a)[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
rz=@(a)[cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1];
rotations={eye(3),rz(pi/3),rx(pi/2),rx(-pi/2),ry(pi/2),ry(-pi/2), ...
    rz(.37)*ry(-.62)*rx(.21),rz(-1.17)*ry(.48)*rx(-.83)};
end
