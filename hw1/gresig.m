function [Msig,Mss] = gresig(flip,T1,T2,TE,TR,dfreq)


Nspin = 200;
M = zeros(3,Nspin);
phi = ([1:Nspin]/Nspin-0.5 ) * 4*pi;


for k=1:Nspin

    %Rotation about y
    Ry = [cos(flip) 0 sin(flip);0 1 0;-sin(flip) 0 cos(flip)];
    [Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
    [Ate,Bte] = freeprecess(TE,T1,T2,dfreq);

    % 	To add the gradient spoiler twist, we just
    %	multiply Atr by zrot(phi):
    Rz = [cos(phi(k)) -sin(phi(k)) 0;sin(phi(k)) cos(phi(k)) 0; 0 0 1];
    Atr = Rz*Atr;
    % Let 	M1 be the magnetization just before the tip.
    %	M2 be just after the tip.
    %	M3 be at TE.
    %
    % then
    %	M2 = Rflip * M1
    %	M3 = Ate * M2 + Bte
    %	M1 = Atr * M3 + Btr
    %
    % Solve for M3...
    %
    %	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)
    
    Mss1 = inv(eye(3)-Ate*Ry*Atr) * (Ate*Ry*Btr+Bte);
    Msig1 = Mss1(1)+i*Mss1(2);
   M(:,k)=Mss1;
end


Mss = mean(M')';
Msig = Mss(1)+i*Mss(2);
