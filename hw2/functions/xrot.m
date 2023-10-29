function Rx=xrot(phi)
% function Rx=xrot(phi)
% rotate magetizatoin vector about x by angle phi

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];