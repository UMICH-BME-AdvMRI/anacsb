%	function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)
%
%	Function calculates the signal from an RF-spoiled sequence
%	after Nex excitations.
%
function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)

if (nargin < 8)
	inc = 117/180*pi;
end;
if (nargin < 7)
	Nex = 100;
end;
if (nargin < 6)
	df = 0;
end;

Nf = 100;	% Simulate 100 different gradient-spoiled spins.
phi = [1:Nf]/Nf*2*pi;

M=zeros(3,Nf,Nex+1);

	
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);

M = [zeros(2,Nf);ones(1,Nf)];
on = ones(1,Nf);
	
Rfph = 0;
Rfinc = inc;

for n=1:Nex

    Rz = [cos(flip) -sin(flip) 0;sin(flip) cos(flip) 0; 0 0 1];
    Rx = [1 0 0; 0 cos(flip) -sin(flip);0 sin(flip) cos(flip)];
    Rth = inv(Rz)*Rx*Rz;
	A = Ate * Rth;
	B = Bte;
	M = A*M+B*on;

	Msig = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph);
	Mss = M;

	M=Atr*M+Btr*on;

	for k=1:Nf
        Rz = [cos(phi(k)) -sin(phi(k)) 0;sin(phi(k)) cos(phi(k)) 0; 0 0 1];
		M(:,k) = Rz*M(:,k);
	end;

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end;