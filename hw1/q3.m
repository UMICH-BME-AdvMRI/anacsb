%% Slice Profile simulation
%% 1)
close all;clear all;clc;
TBW = 8;                % time bandwidth product
S = 180;                % mT/(m*ms) ---- mt/(m.s)
tau_rf = 2;             % ms
gammabar = 42.58;       % MHz/T
BW = TBW/tau_rf;        % kHz
delta_z = 5e-3;         % m

Gss = BW/(gammabar*1e3*delta_z);   % T/m
tau_ss = (Gss/S)*1e3;         % ms
total_gss = (2*tau_ss) + tau_rf;   % ms

Apos = (tau_ss*Gss) + (tau_rf*Gss); % T/m ms
Aneg = Apos/2;
tau_ss2 = 25/S;         % ms
tau_rph = (Aneg - (25*1e-3*tau_ss2))/(25*1e-3);

fprintf('Amplitude of the slice-select gradient: %.3f mT/m\n',Gss*1e3);
fprintf('Slice-select gradient total time: %.3f ms\n',total_gss);
%% 2
T1 = 1000;                          % ms
T2 = 100;                           % ms
dt = 0.01;                          % ms, step size
N = 3/dt;                           % total num of sample for simulation
nn = tau_rf/dt;                     % number of samples for RF pulse    
sincperiod = tau_rf/4;               % sinc stretch parameter
timerf  =  [-(nn-1)/2:(nn-1)/2]*dt;   % ms
timetotal = [0:N-1]*dt;
rf = hanning(length(timerf))'.*sinc(timerf./sincperiod); % sinc() is truncated
% rf amplitude, 90 degree pulse: b1
rf90 = (pi/2)/(2*pi*gammabar*1e3*sum(rf)*tau_rf);   % T
rf_signal = zeros(1,N);
rf_signal(1:nn) = rf90*rf;                          % T 
% gss = fftshift(fft(rf_signal));
% gradient slice selection
gss = zeros(1,N);
gss(1,1:nn) = [Gss*ones(1,nn)];        % T/m
% gss(nn+1:nn+1+floor(tau_rph/dt)) = -.025;

pos = -10:dt:10;                   %mm
df = [0 200];                       % Off- resonance frequencies

for jj = 1: length(df)
    M =[0;0;1];	                    % Starting magnetization.
    [A,B] = freeprecess(tau_rf/2,T1,T2,df(jj));
    for ii=1:length(pos)
        for k = 1:length(rf)
	        M = A*M+B;
	        gssrot = 2*pi*gammabar*(pos(ii))*gss(k)*dt/2;
            Rz = [cos(gssrot) -sin(gssrot) 0;sin(gssrot) cos(gssrot) 0; 0 0 1];
	        M = Rz*M;
        
            M = throt(abs(rf_signal(k)),angle(rf_signal(k))) * M;	% RF Rotation.
        
	        M = A*M+B;
	        gssrot = 2*pi*gammabar*(pos(ii))*gss(k)*dt/2;
            Rz = [cos(gssrot) -sin(gssrot) 0;sin(gssrot) cos(gssrot) 0; 0 0 1];
    
	        M = Rz*M;
        end
    
        m(:,ii,jj) = M;
        msig(jj,ii) = M(1)+i*M(2);
    end
end

figure,
subplot(321);
stem(timetotal,rf_signal*1e6);xlabel('Time (ms)'); ylabel('(uT)');title('RF waveform');grid on;
subplot(322);
plot(timetotal,gss);xlabel('Time (ms)'); ylabel('Gss (mT/m)');title('Slice selective gradient');grid on;

subplot(323);
plot(pos,squeeze(m(1,:,1)));xlabel('Position (mm)');ylabel('M_{x}');title('M_{x}');hold on;grid on;
plot(pos,squeeze(m(1,:,2)));xlabel('Position (mm)');legend('0 Hz off-ressonance','200 Hz off-resonnace');

subplot(324);
plot(pos,squeeze(m(2,:,1)));xlabel('Position (mm)');ylabel('M_{y}');title('M_{y}');hold on;grid on;
plot(pos,squeeze(m(2,:,2)));legend('0 Hz off-ressonance','200 Hz off-resonnace');

subplot(325);
plot(pos,squeeze(m(3,:,1)));xlabel('Position (mm)');ylabel('M_{z}');title('M_{z}');hold on;grid on;
plot(pos,squeeze(m(3,:,2)));legend('0 Hz off-ressonance','200 Hz off-resonnace');

subplot(326);
plot(pos,abs(msig(1,:)));xlabel('Position (mm)');ylabel('M_{xy}');title('M_{xy}');hold on;grid on;
plot(pos,abs(msig(2,:)));legend('0 Hz off-ressonance','200 Hz off-resonnace'); set(gcf,'color','w');
fontsize(14,"points");
%% 3
% rf amplitude, 30 degree pulse: b1
rf30 = (pi/6)/(2*pi*gammabar*1e3*sum(rf)*tau_rf);   % T
rf30_signal = zeros(1,N);
rf30_signal(1:nn) = rf30*rf;   
pos = -10:dt:10;                  %mm
df = [0 200];
[msig30,m30] = sprofile(rf30_signal,gss,tau_rf,T1,T2,pos,df,dt,gammabar);

figure,
subplot(321);
stem(timetotal,rf30_signal*1e6);xlabel('Time (ms)'); ylabel('(uT)');title('RF waveform 30°');grid on;
subplot(322);
plot(timetotal,gss);xlabel('Time (ms)'); ylabel('Gss (mT/m)');grid on;

subplot(323);
plot(pos,squeeze(m30(1,:,1)));xlabel('Position (mm)');ylabel('M_{x}');title('M_{x}');hold on;
plot(pos,squeeze(m30(1,:,2)));xlabel('Position (mm)');legend('0 Hz off-ressonance','200 Hz off-resonnace');grid on;

subplot(324);
plot(pos,squeeze(m30(2,:,1)));xlabel('Position (mm)');ylabel('M_{y}');title('M_{y}');hold on;grid on;
plot(pos,squeeze(m30(2,:,2)));legend('0 Hz off-ressonance','200 Hz off-resonnace');grid on;

subplot(325);
plot(pos,squeeze(m30(3,:,1)));xlabel('Position (mm)');ylabel('M_{z}');title('M_{z}');hold on
plot(pos,squeeze(m30(3,:,2)));legend('0 Hz off-ressonance','200 Hz off-resonnace');grid on;

subplot(326);
plot(pos,abs(msig30(1,:)));xlabel('Position (mm)');ylabel('M_{xy}');title('M_{xy}');hold on
plot(pos,abs(msig30(2,:)));legend('0 Hz off-ressonance','200 Hz off-resonnace');grid on;
fontsize(14,"points");

% RF 10 degree pulse
rf10 = (pi/18)/(2*pi*gammabar*1e3*sum(rf)*tau_rf);   % T
rf10_signal = zeros(1,N);
rf10_signal(1:nn) = rf10*rf;   
pos = -10:dt:10;                   %mm
df = 0;
[msig10,m10] = sprofile(rf10_signal,gss,tau_rf,T1,T2,pos,df,dt,gammabar);

figure,
subplot(321);
stem(timetotal,rf10_signal*1e6);xlabel('Time (ms)'); ylabel('(uT)');title('RF waveform 10°');grid on;
subplot(322);
plot(timetotal,gss);xlabel('Time (ms)'); ylabel('Gss (mT/m)');grid on;

subplot(323);
plot(pos,squeeze(m10(1,:)));xlabel('Position (mm)');ylabel('M_{x}');title('M_{x}');hold on;grid on;

subplot(324);
plot(pos,squeeze(m10(2,:)));xlabel('Position (mm)');ylabel('M_{y}');title('M_{y}');grid on;

subplot(325);
plot(pos,squeeze(m10(3,:)));xlabel('Position (mm)');ylabel('M_{z}');title('M_{z}');grid on;

subplot(326);
plot(pos,abs(msig10(1,:)));xlabel('Position (mm)');ylabel('M_{xy}');title('M_{xy}');grid on;
fontsize(14,"points");

% Comparison FFT
figure,
subplot(221);
plot(pos,abs(msig(1,:)));xlabel('Position (mm)');ylabel('M_{xy}');title('M_{xy}, 90°');
subplot(222);
plot(pos,abs(msig10));title('M_{xy}, 10°');grid on;xlabel('Position (mm)');ylabel('M_{xy}');
% with FFT
subplot(223);
plot(pos,abs(ifftshift(fft(ifftshift(msig(1,:))))));title('FFT M_{xy}, 90°');grid on;
xlabel('Position (mm)');ylabel('M_{xy}');
subplot(224);
plot(pos,abs(ifftshift(fft(ifftshift(msig10)))));title('FFT M_{xy}, 10°');grid on;
xlabel('Position (mm)');ylabel('M_{xy}');fontsize(14,"points");

% Changing T2
T2 = 2;             % ms
[msig90,m90] = sprofile(rf_signal,gss,tau_rf,T1,T2,pos,df,dt,gammabar);
[msig30,m30] = sprofile(rf30_signal,gss,tau_rf,T1,T2,pos,df,dt,gammabar);
[msig10,m10] = sprofile(rf10_signal,gss,tau_rf,T1,T2,pos,df,dt,gammabar);

figure,
subplot(321);
plot(timetotal,rf_signal*1e6);xlabel('Time (ms)'); ylabel('(uT)');title('RF waveform');hold on;
plot(timetotal,rf30_signal*1e6);grid on;
plot(timetotal,rf10_signal*1e6); legend('90°','30°','10°');

subplot(322);
plot(timetotal,gss);xlabel('Time (ms)'); ylabel('Gss (mT/m)');grid on;

subplot(323);
plot(pos,squeeze(m90(1,:)));xlabel('Position (mm)');ylabel('M_{x}');title('M_{x}');hold on
plot(pos,squeeze(m30(1,:)));grid on;
plot(pos,squeeze(m10(1,:)));legend('90°','30°','10°');

subplot(324);
plot(pos,squeeze(m90(2,:)));xlabel('Position (mm)');ylabel('M_{y}');title('M_{y}');hold on
plot(pos,squeeze(m30(2,:)));grid on;
plot(pos,squeeze(m10(2,:)));legend('90°','30°','10°');

subplot(325);
plot(pos,squeeze(m90(3,:)));xlabel('Position (mm)');ylabel('M_{z}');title('M_{z}');hold on
plot(pos,squeeze(m30(3,:)));grid on;
plot(pos,squeeze(m10(3,:)));legend('90°','30°','10°');

subplot(326);
plot(pos,abs(msig90(1,:)));xlabel('Position (mm)');ylabel('M_{xy}');title('M_{xy}');hold on
plot(pos,abs(msig30(1,:)));grid on;
plot(pos,abs(msig10(1,:)));legend('90°','30°','10°');fontsize(14,"points");
%% 4
T1 = 1000;                          % ms
T2 = 100;                           % ms

pos = -10:dt:10;                    %mm
df = 0;                             % Off- resonance frequencies

for ii=1:length(pos)
    M =[0;0;1];	                    % Starting magnetization.
    M = throt(abs(rf_signal(1)),angle(rf_signal(1))) * M;	% RF Rotation.

    for k = 2:length(rf)
            [A,B] = freeprecess(dt,T1,T2,0);
	        M = A*M+B;
	        gssrot = 2*pi*gammabar*(pos(ii))*gss(k)*dt/2;
            Rz = [cos(gssrot) -sin(gssrot) 0;sin(gssrot) cos(gssrot) 0; 0 0 1];
	        M = Rz*M;
        
            M = throt(abs(rf_signal(k)),angle(rf_signal(k))) * M;	% RF Rotation.
     end
    
        m(:,ii) = M;
        msig(ii) = M(1)+i*M(2);
end

figure,
subplot(211);plot(pos,abs(msig(1,:)));xlabel('position (mm)');ylabel('Amplitude');
subplot(212);plot(pos,angle(msig(1,:)));xlabel('position (mm)');ylabel('Phase(rad)');

%% 5
T1 = 1000;                          % ms
T2 = 100;                           % ms
dt = 0.01;                          % ms, step size
N = 3/dt;                           % total num of sample for simulation
nn = tau_rf/dt;                     % number of samples for RF pulse    
sincperiod = tau_rf/4;               % sinc stretch parameter
timerf  =  [-(nn-1)/2:(nn-1)/2]*dt;   % ms
timetotal = [0:N-1]*dt;
rf = hanning(length(timerf))'.*sinc(timerf./sincperiod); % sinc() is truncated
[~,ind]= find(rf==max(rf(:)));
gss = zeros(1,N);
gss(1,1:nn) = [Gss*ones(1,nn)];        % T/m
gss(nn+1:nn+1+floor(tau_rph/dt)) = -.025;         % T/m

dz = 0.005;                              % m

rfsig = rf;
timerf2  = [0:nn-1]*dt;   % ms
% rfmod = exp(i*2*pi*gammabar*1e3*gss(2)*dz.*timerf2);
% rf2 = [rfsig.*rfmod + rfsig.*conj(rfmod)];
% rfsig = rf2;

rf3 = 0;
% adding the phase 
phi = 0:pi/2:2*pi;                      % rad
z = [-40:20:40];                        % mm      

P=0;
for kk = 1: length(z)
    rfmodsum = exp(i*2*pi*gammabar*gss(2)*z(kk).*timerf2);
    phse = phase(rfmodsum(ind(1)));
    ptmp = rfmodsum.*exp(phse);
    % rf2 = [rfsig.*rfmod + rfsig.*conj(rfmod)]; 
    % rfamp = phi(kk)/(2*pi*gammabar*1e3*sum(rf2)*tau_rf);
    P = P + ptmp;%.*rfamp;
end

rfmb = rfsig.*P;
rf90 = (pi/2)/(2*pi*gammabar*1e3*sum(rfmb)*tau_rf);   % T

% rf amplitude, 90 degree pulse: b1
% rf90 = (pi/2)/(2*pi*gammabar*1e3*sum(rfsig)*tau_rf);   % T

rf_signal = zeros(1,N);
% rf_signal(1:nn) = rf2.*rf90;                          % T 

rf_signal(1:nn) = rfmb.*rf90;                          % T 


df = 0;
pos = [-50:dt:50];                              % mm

[mmsig,mm] = sprofile(rf_signal,gss,tau_rf,T1,T2,pos,df,dt,gammabar);

figure;
plot(timetotal,abs(rf_signal)*1e6);xlabel('Time (ms)'); 
ylabel('(uT)');title('SMS pulse');grid on;fontsize(14,"points");
figure;
plot(pos,abs(mmsig));title('Slice Profile');grid on;
xlabel('Position (mm)');ylabel('Signal Magnitude');grid on;fontsize(14,"points");

