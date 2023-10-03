%% Balanced and Spoiled Steady-State sequences
%% 2a)
% Reference: http://mrsrl.stanford.edu/~brian/bloch/
clear all;clc;close all;
T1 = 1000;	% ms.
T2 = 100;	% ms.
TR = [5 10 20];	% ms.
TE = TR/2;	% ms.

flip = pi/3;    %radians
df = [-200:200]; 	%Hz
S = zeros(length(df),length(TE));

for n=1:length(TE)
	for k=1:length(df)
        %rotation in y:
        Rflip = [cos(flip) 0 sin(flip);0 1 0;-sin(flip) 0 cos(flip)];
        
        [Atr,Btr] = freeprecess(TR(n)-TE(n),T1,T2,df(k));
        [Ate,Bte] = freeprecess(TE(n),T1,T2,df(k));
        % at TR:
        Mtr = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
		S(k,n)=Mtr(1)+i*Mtr(2);
    end
end
% Plots
figure,plot(df,abs(S));
xlabel('Frequency (Hz)');ylabel('Magnitude');
title('Frequency response of a bSSFP');
grid on; set(gcf,'color','w');
legend('TR = 5', 'TR = 10', 'TR = 20');
fontsize(14,"points")
% amplitude agrees with the periodic resonant freq (1/TR) and it is flat
% when TE/2
%% 2b) i
clear all;clf;
T1 = 1000;	% ms.
T2 = 100;	% ms.
TR = 10;	% ms.
TE = TR/2;	% ms.
flip = (10/180)*pi;    %radians

df = [-100:5:100]; 	%Hz

Sig = zeros(length(df),length(TE));

for n=1:length(TE)
	for k=1:length(df)
		[Msig,Mss] = gresig(flip,T1,T2,TE(n),TR(n),df(k));
		Sig(k,n)=Mss(1)+i*Mss(2);
    end

end
%plot the Results 
figure;plot(df,abs(Sig));
xlabel('Frequency (Hz)');ylabel('|M_{ss}|');
grid on;legend('TR = 10');ylim([0 .1]);
title('M_{ss} vs. Frequency');set(gcf,'color','w');

%% 2b) ii
clear all;clf;

T1 = 1000;	% ms.
T2 = 100;	% ms.
TR = 10;	% ms.
TE = 5;	% ms.

flip = (10/180)*pi;    %radians
inc = pi;

num_exc = 100;
num_spoil = 100;
% num_exc = 4;
% num_spoil = 4;
phi = [1:num_spoil]/num_spoil*2*pi;
% phi = 2*[pi 2*pi 4*pi 8*pi];
M=zeros(3,num_spoil,num_exc+1);
Msig = zeros(1,num_exc);

[Ate,Bte] = freeprecess(TE,T1,T2,0);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,0);

M = [zeros(2,num_spoil);ones(1,num_spoil)];
on = ones(1,num_spoil);
	
Rfph = 0;
Rfinc = inc;

for n=1:num_exc

    Rz = [cos(flip) -sin(flip) 0;sin(flip) cos(flip) 0; 0 0 1];
    Rx = [1 0 0; 0 cos(flip) -sin(flip);0 sin(flip) cos(flip)];
    Rth = inv(Rz)*Rx*Rz;

	A = Ate * Rth;
	B = Bte;
	M = A*M+B*on;

	Msig(n) = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph);
	
	M=Atr*M+Btr*on;

	for k=1:num_spoil
        Rz = [cos(phi(k)) -sin(phi(k)) 0;sin(phi(k)) cos(phi(k)) 0; 0 0 1];
		M(:,k) = Rz*M(:,k);
    end

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end

%plot
time = [0:num_exc-1]*TR+TE;
plot(phi,abs(Msig));xlabel('dephasing moment (rad)'); title('TR=10 ms, TE = 5 ms');
ylabel('Magnitude');grid on;
%% 2b) iii
clear all; clf;
df = 0;		% Hz off-resonance.
T1 = 1000;	% ms.
T2 = 100;	% ms.
TE = 5;		% ms.
TR = 10;
flip = [0:0.01:1]*pi;	% radians.
inc = 117/180*pi; % provides sufficient large phase increment 
Nex = 100;

for k=1:length(flip)
	sig1(k)=spgrsignal(flip(k),T1,T2,TE,TR,df,Nex,inc);
	[Msig,M]=srsignal(flip(k),T1,T2,TE,TR,df);
	sig2(k)=M(1)+i*M(2);
end
figure;
plot(flip,abs(sig1));hold on;plot(flip,abs(sig2));
xlabel('Flip (rad)');ylabel('Magnitude');
grid on;title('Gradient spoiler and RF spoling signals');
legend('RFspoiling','Gradient spoiler');

% when we simulate both effects, we need to find the RF phase that best
% eliminates the transverse magnetization, which corresponds to the minimum
% steady state signal
% Find optimal phase
[min_signal, min_idx] = find(sig1(1,3:end) == min(sig1(1,3:end)));
optimal_phase = flip(min_idx);

fprintf('Optimal RF phase: %.3f \n',optimal_phase);
