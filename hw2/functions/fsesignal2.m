
function [Msig,Mss] = fsesignal2(T1,T2,TE,TR,dfreq,ETL)

Rflip = yrot(pi/2);	% Rotation from Excitation  (usually 90)
Rrefoc = xrot(pi);	% Rotation from Refocusing (usually 180)

[Atr,Btr] = freeprecess(TR-ETL*TE,T1,T2,dfreq);	% Propagation last echo to TR
[Ate2,Bte2] = freeprecess(TE/2,T1,T2,dfreq);	% Propagation over TE/2

% Neglect residual transverse magnetization prior to excitation.
Atr = [0 0 0;0 0 0;0 0 1]*Atr;	% Retain only Mz component.


% Since ETL varies, let's keep a "running" A and B.  We'll
% calculate the steady-state signal just after the tip, Rflip.

% Initial.
A=eye(3);
B=[0 0 0]';


% For each echo, we "propagate" A and B:
for k=1:ETL
	A = Ate2*Rrefoc*Ate2*A;			% TE/2 -- Refoc -- TE/2
	B = Bte2+Ate2*Rrefoc*(Bte2+Ate2*B);
end;


% Propagate A and B through to just after flip, and calculate steady-state.
A = Rflip*Atr*A;
B = Rflip*(Btr+Atr*B);

Mss = pinv(eye(3)-A)*B;	% -- Steady state is right after 90 pulse!
M = Mss;


% Calculate signal on each echo.
for k=1:ETL
	M = Ate2*Rrefoc*Ate2*M + Bte2+ Ate2*Rrefoc*Bte2;
	Mss(:,k)=M;
	Msig(k)=M(1)+i*M(2);
end;

