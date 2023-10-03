%% Spatial encoding
%% 1a)
TBW = 2;        
tau_rf = 1*1e-3;            % s
delta_z = 3*1e-3;           % m
gamma_bar = 42.58*1e6;      % Hz/T
S = 180;                    % T/m/s

%Calculate Gss
Gss = TBW/(tau_rf*delta_z*gamma_bar);       % T/m

%Calculate tau_ss
tau_ss_rise = Gss/S;                     % s
total_ss = (2*tau_ss_rise) + tau_rf;     % s
% Total area
A2 = (tau_ss_rise*Gss) + (Gss*tau_rf);
A = A2/2;
Grp = A/(total_ss/2);                    % T/m
% Grp = 25*1e-3;                    % T/m
tau_rp_rise = Grp/S;                     % s
tau_rp = total_ss/2;
total_ss = total_ss*1e3;                % ms
fprintf('Amplitude of the slice-select gradient: %.3f mT/m\n',Gss*1e3);
fprintf('Amplitude of the slice-rephasing gradient: %.3f mT/m\n',Grp*1e3);
fprintf('Slice-select gradient total time: %.3f ms\n',total_ss);
%% 1b)
tau_sr = total_ss/2;
N = 256;                   % number of samples
delta_y = 1.2*1e-3;        % m
Gmax_pe = 25*1e-3;         % T/m
% Calculate  tau_pe = 1 / (gamma_bar * Gmax_pe * 2 * delta_y)
tau_pe = 1 / (gamma_bar * Gmax_pe * 2 * delta_y);    % s

%Calculate tau_pe_rise
tau_pe_rise = Gmax_pe/S;                    % s
total_pe = (2*tau_pe_rise) + tau_pe;        % s
fprintf('Phase-encoding gradient total time %.3f ms\n',total_pe*1e3);
%% 1c)
rBW = 750;                  % Hz/pixel
delta_x = 1.2*1e-3;         % m
Ts = 1/rBW;                     % s

%Calculate total rBW
totalRBW = rBW*N;               % Hz

% Calculate FOVx
kx_max = 1/(2*delta_x);         % m^-1
kx = 2*kx_max;                  % m^-1
delta_kx = kx/N;                % m^-1
FOVx = 1/delta_kx;              % m^-1

% Calculate Gread = (totalRBW)/(gamma_bar*FOV)
Gread = (totalRBW)/(gamma_bar*FOVx);    % T/m
% Calculate tau_ro_rise
tau_ro_rise = Gread/S;                  % s
total_ro = (2*tau_ro_rise) + Ts;        % s
% Total area
A_2 = (tau_ro_rise*Gread) + (Gread*Ts);
A = A_2/2;
Gpre_ro = A/(total_ro/2);                    % T/m
tau_pre_ro_rise = Gpre_ro/S;                     % s

% total time
Ttime = total_ss +  1.5*total_ro*1e3;
fprintf('Read-out gradient total time %.3f ms\n',total_ro*1e3);
fprintf('Total time %.3f ms\n',Ttime);
fprintf('TE= %.3f ms\n', Ttime - (total_ro*1e3/2));
fprintf('TR= %.3f ms\n', Ttime)
%% 1d)
dt = 0.01;           % ms
time = 0: dt:Ttime;       % ms
tmp = round(tau_ss_rise*1e3,4);
n = floor(tmp/dt);
npoints = tau_rf*1e3/dt;
t_sinc = tau_rf*1e3/8;    % ms
rftime = [-(npoints-1)/2:(npoints-1)/2]*dt;
rf = hanning(npoints)'.*sinc(rftime/t_sinc);

rfpulse = zeros(1,length(time));
rfpulse(floor(tmp/dt):floor(tmp/dt) + floor(tau_rf*1e3/dt)-1) = rf;

figure,
subplot(411),plot(time,rfpulse);ylabel('RF pulse');grid on;

subplot(412);ylabel('G_z (mT/m)');hold on;
plot([0 tau_ss_rise*1e3],[0 Gss*1e3],"red");hold on;
plot([tau_ss_rise*1e3,tau_ss_rise*1e3 + tau_rf*1e3],[Gss*1e3,Gss*1e3],'r');
plot([tau_ss_rise*1e3 + tau_rf*1e3,total_ss ],[Gss*1e3 0],'r');
plot([total_ss, total_ss+tau_rp_rise*1e3], [0 -Grp*1e3],'r');
plot([total_ss+tau_rp_rise*1e3, 1.5*total_ss-tau_rp_rise*1e3],[-Grp*1e3 -Grp*1e3],'r');
plot([1.5*total_ss-tau_rp_rise*1e3, 1.5*total_ss],[-Grp*1e3 0],'r');
plot([1.5*total_ss, Ttime],[0 0],'r');grid on;

subplot(413);
plot([0 total_ss],[0 0],'b');hold on;
plot([total_ss (total_ss+ tau_pe_rise*1e3)],[0 Gmax_pe*1e3],'b');
plot([(total_ss+ tau_pe_rise*1e3),(total_ss+ tau_pe_rise*1e3)+(tau_pe*1e3)],[Gmax_pe*1e3,Gmax_pe*1e3],'b');
plot([(total_ss+ tau_pe_rise*1e3),(total_ss+ tau_pe_rise*1e3)+(tau_pe*1e3)],[Gmax_pe*1e3,Gmax_pe*1e3],'b');
plot([(total_ss+ tau_pe_rise*1e3)+(tau_pe*1e3),total_ss + total_pe*1e3],[Gmax_pe*1e3 0],'b')
plot([total_ss + total_pe*1e3, Ttime],[0 0],'b');ylabel('G_y (mT/m)');grid on;

subplot(414);
plot([0,(total_ss)],[0 0],'k');hold on;
plot([total_ss, total_ss+tau_pre_ro_rise*1e3],[0 -Gpre_ro*1e3],'k');
plot([total_ss+tau_pre_ro_rise*1e3,total_ss+(total_ro/2)*1e3- tau_pre_ro_rise*1e3],[-Gpre_ro*1e3 -Gpre_ro*1e3],'k');
plot([total_ss+(total_ro/2)*1e3- tau_pre_ro_rise*1e3,total_ss+(total_ro/2)*1e3 ],[-Gpre_ro*1e3 0],'k');
plot([(total_ss + (total_ro/2)*1e3), total_ss+(total_ro/2)*1e3+tau_ro_rise*1e3],[0 Gread*1e3],'k');
plot([total_ss+(total_ro/2)*1e3+tau_ro_rise*1e3 , total_ss+(total_ro/2)*1e3+tau_ro_rise*1e3+Ts*1e3],[Gread*1e3 Gread*1e3],'k');hold on;
plot([total_ss+(total_ro/2)*1e3+tau_ro_rise*1e3+Ts*1e3,total_ss+(total_ro/2)*1e3+total_ro*1e3],[Gread*1e3,0],'k')
plot([total_ss+(total_ro/2)*1e3+total_ro*1e3, Ttime],[0 0],'k');ylabel('G_x (mT/m)');xlabel('Time (ms)');
grid on;
fontsize(14,"points")

fprintf('Ways to decrease TR and TE:\n');
fprintf('(1) Acquire half of the k-space\n');
fprintf('(2) Smaller flip angle\n');
fprintf('(3) Use parallel imaging\n');
