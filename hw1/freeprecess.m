function [Afp,Bfp]=freeprecess(T,T1,T2,df)
%   function [Afp,Bfp]=freeprecess(T,T1,T2,df)
%	Function simulates free precession and decay
%	over a time interval T, given relaxation times T1 and T2 (ms)
%	and off-resonance df.  Times in ms, off-resonance (Hz)

phi = 2*pi*df*T*10^-3;	% Resonant precession, radians.
E1 = exp(-T/T1);        % T1 relaxation	
E2 = exp(-T/T2);        % T2 relaxation

%Calculate rotation around z:
Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];

% M = A*Rz + B
Afp = [E2 0 0;0 E2 0;0 0 E1]*Rz;
Bfp = [0 0 1-E1]';
