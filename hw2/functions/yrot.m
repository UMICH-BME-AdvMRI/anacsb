function Ry=yrot(phi)
% function Ry=yrot(phi)
% rotate magetizatoin vector about y by angle phi

Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];