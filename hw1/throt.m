function Rth=throt(phi,theta)
% function Rth=throt(phi,theta)
% rotate magnetization vector by some angle theta

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
Rz = [cos(-theta) -sin(-theta) 0;sin(-theta) cos(-theta) 0; 0 0 1];

Rth = inv(Rz)*Rx*Rz;

end
