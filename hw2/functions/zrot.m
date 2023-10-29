function Rz=zrot(phi)
% function Rz=zrot(phi)
% rotate magetizatoin vector about z by angle phi
% assume right-handed coordinates and rotations

Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];
