%% Problem 1: EPG
addpath(('functions'));
%% a)
clear all;clc;
disp('-- Spin Echo with 180y pulses --')
etl = 64;
esp = 5*1e-3;    % sec
TR = etl*esp;
N = 2*etl;
Qall = zeros(4*etl,etl);    % store f, F* states
Zall = zeros(2*etl,etl);    % store Z states
T1 = [200:100:1500]*1e-3;              % sec
T2 = [50:30:300]*1e-3;                 % sec
timet = [1:etl]*esp*1000;

numcomb = length(T1)*length(T2);
signal = zeros(numcomb,etl);
ff = 1;
alpha = 180;
flipangle = pi/180*[alpha+(90-alpha/2) alpha*ones(1,etl-1)];
for ii = 1 :  length(T1)
    for jj = 1 : length(T2)
        Q = epg_m0(N); fprintf('T1 = %.2f ms, T2 = %.2f ms \n',T1(ii)*1e3,T2(jj)*1e3);
        Q = epg_relax(Q,T1(ii),T2(jj),esp); 
        % 90 RF pulse
        Q = epg_rf(Q,pi/2,pi/2); 
            for ech = 1 : etl
                % Refoc. pulses
                Q = epg_grelax(Q,T1(ii),T2(jj),esp/2,1,0,1,1);  % -- Left crusher
                Q = epg_rf(Q,abs(flipangle(ech)),angle(flipangle(ech))); % -- Refoc. RF
                Q = epg_grelax(Q,T1(ii),T2(jj),esp/2,1,0,1,1);   % -- Right crusher
            
                signal(ff,ech) = Q(1,1);
                Qall(2*etl:4*etl-1,ech) = Q(2,:).';	        %negative states
                Qall(1:2*etl,ech) = flipud(Q(1,:).');  % Put in positive, overwrite center.
                Zall(:,ech) = Q(3,:).';
            end
        ff = ff + 1;
    end
end

%% a.i.)
selT1T2 = [1 75 67 123 126];
ii = 1 :  length(T1);
jj = 1 : length(T2);
C = combinations(ii,jj);
figure('Renderer', 'painters', 'Position', [50 50 900 600]);set(gcf,'DefaultLineLineWidth',2);
plot(timet,abs(signal(selT1T2(1),:))); hold on;
plot(timet,abs(signal(selT1T2(2),:))); hold on;
plot(timet,abs(signal(selT1T2(3),:))); hold on;
plot(timet,abs(signal(selT1T2(4),:))); hold on;
plot(timet,abs(signal(selT1T2(5),:))); 
legend1 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(1)))*1e3,T2(C.jj(selT1T2(1)))*1e3);
legend2 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(2)))*1e3,T2(C.jj(selT1T2(2)))*1e3);
legend3 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(3)))*1e3,T2(C.jj(selT1T2(3)))*1e3);
legend4 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(4)))*1e3,T2(C.jj(selT1T2(4)))*1e3);
legend5 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(5)))*1e3,T2(C.jj(selT1T2(5)))*1e3);
fontsize(gcf, 24, "points");
legend(legend1, legend2, legend3,legend4,legend5);
lplot('Echo Time (ms)','Signal','Signal vs. Echo Time, \alpha = 180');
if ~exist('/figs_q1', 'dir')
       mkdir('figs_q1')
end
saveas(gcf,'./figs_q1/a_i.png');
close;

%% a.ii)
clear all;clc;close all;
disp('-- Spin Echo with 120y pulses --')
etl = 64;
esp = 5*1e-3;    % sec
TR = etl*esp;
N = 2*etl;
Qall = zeros(4*etl,etl);    % store f, F* states
Zall = zeros(2*etl,etl);    % store Z states
T1 = [200:100:1500]*1e-3;              % sec
T2 = [50:30:300]*1e-3;                 % sec
timet = [1:etl]*esp*1000;

numcomb = length(T1)*length(T2);
signal = zeros(numcomb,etl);
ff = 1;
alpha = 120;
flipangle = pi/180*[alpha+(90-alpha/2) alpha*ones(1,etl-1)];

for ii = 1 :  length(T1)
    for jj = 1 : length(T2)
        Q = epg_m0(N); fprintf('T1 = %.2f ms, T2 = %.2f ms \n',T1(ii)*1e3,T2(jj)*1e3);
        Q = epg_relax(Q,T1(ii),T2(jj),esp); 
        % 90 RF pulse
        Q = epg_rf(Q,pi/2,pi/2); 
            for ech = 1 : etl
                % Refoc. pulses
                Q = epg_grelax(Q,T1(ii),T2(jj),esp/2,1,0,1,1);  % -- Left crusher
                Q = epg_rf(Q,abs(flipangle(ech)),angle(flipangle(ech))); % -- Refoc. RF
                Q = epg_grelax(Q,T1(ii),T2(jj),esp/2,1,0,1,1);   % -- Right crusher
            
                signal(ff,ech) = Q(1,1);
                Qall(2*etl:4*etl-1,ech) = Q(2,:).';	        %negative states
                Qall(1:2*etl,ech) = flipud(Q(1,:).');  % Put in positive, overwrite center.
                Zall(:,ech) = Q(3,:).';
            end
        ff = ff + 1;
    end
end

selT1T2 = [1 75 67 123 126];
ii = 1 :  length(T1);
jj = 1 : length(T2);
C = combinations(ii,jj);
figure('Renderer', 'painters', 'Position', [50 50 900 600]);set(gcf,'DefaultLineLineWidth',2);
plot(timet,abs(signal(selT1T2(1),:))); hold on;
plot(timet,abs(signal(selT1T2(2),:))); hold on;
plot(timet,abs(signal(selT1T2(3),:))); hold on;
plot(timet,abs(signal(selT1T2(4),:))); hold on;
plot(timet,abs(signal(selT1T2(5),:))); 
legend1 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(1)))*1e3,T2(C.jj(selT1T2(1)))*1e3);
legend2 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(2)))*1e3,T2(C.jj(selT1T2(2)))*1e3);
legend3 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(3)))*1e3,T2(C.jj(selT1T2(3)))*1e3);
legend4 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(4)))*1e3,T2(C.jj(selT1T2(4)))*1e3);
legend5 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(5)))*1e3,T2(C.jj(selT1T2(5)))*1e3);
fontsize(gcf, 24, "points");
legend(legend1, legend2, legend3,legend4,legend5);
lplot('Echo Time (ms)','Signal','Signal vs. Echo Time, \alpha = 120');
if ~exist('/figs_q1', 'dir')
       mkdir('figs_q1')
end
saveas(gcf,'./figs_q1/a_ii.png');
close;

%% a.iii)
clear all;clc;
disp('-- Spin Echo with 60y pulses --')
etl = 64;
esp = 5*1e-3;    % sec
TR = etl*esp;
N = 2*etl;
Qall = zeros(4*etl,etl);    % store f, F* states
Zall = zeros(2*etl,etl);    % store Z states
T1 = [200:100:1500]*1e-3;              % sec
T2 = [50:30:300]*1e-3;                 % sec
timet = [1:etl]*esp*1000;

numcomb = length(T1)*length(T2);
signal = zeros(numcomb,etl);
ff = 1;
alpha = 60;
flipangle = pi/180*[alpha+(90-alpha/2) alpha*ones(1,etl-1)];

for ii = 1 :  length(T1)
    for jj = 1 : length(T2)
        Q = epg_m0(N); fprintf('T1 = %.2f ms, T2 = %.2f ms \n',T1(ii)*1e3,T2(jj)*1e3);
        Q = epg_relax(Q,T1(ii),T2(jj),esp); 
        % 90 RF pulse
        Q = epg_rf(Q,pi/2,pi/2); 
            for ech = 1 : etl
                % Refoc. pulses
                Q = epg_grelax(Q,T1(ii),T2(jj),esp/2,1,0,1,1);  % -- Left crusher
                Q = epg_rf(Q,abs(flipangle(ech)),angle(flipangle(ech))); % -- Refoc. RF
                Q = epg_grelax(Q,T1(ii),T2(jj),esp/2,1,0,1,1);   % -- Right crusher
            
                signal(ff,ech) = Q(1,1);
                Qall(2*etl:4*etl-1,ech) = Q(2,:).';	        %negative states
                Qall(1:2*etl,ech) = flipud(Q(1,:).');  % Put in positive, overwrite center.
                Zall(:,ech) = Q(3,:).';
            end
        ff = ff + 1;
    end
end

selT1T2 = [1 75 67 123 126];
ii = 1 :  length(T1);
jj = 1 : length(T2);
C = combinations(ii,jj);
figure('Renderer', 'painters', 'Position', [50 50 900 600]);set(gcf,'DefaultLineLineWidth',2);
plot(timet,abs(signal(selT1T2(1),:))); hold on;
plot(timet,abs(signal(selT1T2(2),:))); hold on;
plot(timet,abs(signal(selT1T2(3),:))); hold on;
plot(timet,abs(signal(selT1T2(4),:))); hold on;
plot(timet,abs(signal(selT1T2(5),:))); 
legend1 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(1)))*1e3,T2(C.jj(selT1T2(1)))*1e3);
legend2 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(2)))*1e3,T2(C.jj(selT1T2(2)))*1e3);
legend3 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(3)))*1e3,T2(C.jj(selT1T2(3)))*1e3);
legend4 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(4)))*1e3,T2(C.jj(selT1T2(4)))*1e3);
legend5 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(C.ii(selT1T2(5)))*1e3,T2(C.jj(selT1T2(5)))*1e3);
fontsize(gcf, 24, "points");
legend(legend1, legend2, legend3,legend4,legend5);
lplot('Echo Time (ms)','Signal','Signal vs. Echo Time, \alpha = 60');
saveas(gcf,'./figs_q1/a_iii.png');
close;
%% b)
clear all;close all;clc;
disp('-- Spin Echo with alpha=180 --')
esp = 5*1e-3;    % sec
T1 = [200:100:1500]*1e-3;              % sec
T2 = [50:30:300]*1e-3;                 % sec

alpha = 180;
num_flips = 64;
flips = pi/180*[alpha+(90-alpha/2) alpha*ones(1,num_flips-1)];
for ii=1:length(T2)
    for jj=1:length(T1)
        S(ii,jj,:) = abs(epg_cpmg(flips,[],T1(jj),T2(ii),esp));
    end
end
sels = [6,16,32,48];

figure('Renderer', 'painters', 'Position', [50 50 1600 1100]);set(gcf,'DefaultLineLineWidth',2);
for ii=1:length(sels)
    subplot(1,4,ii)
    contour(T1*1e3,T2*1e3,S(:,:,sels(ii)),'LineWidth',2)
    title(['Echo #' num2str(sels(ii))])
    xlabel('T1');ylabel('T2');
end
sgtitle('\alpha = 180°') 
fontsize(gcf, 24, "points");
saveas(gcf,'./figs_q1/b_180.png'); close;

disp('-- Spin Echo with alpha=120 --')

alpha = 120;
num_flips = nTR*etl;
flips = pi/180*[alpha+(90-alpha/2) alpha*ones(1,num_flips-1)];
for ii=1:length(T2)
    for jj=1:length(T1)
        S(ii,jj,:) = abs(epg_cpmg(flips,[],T1(jj),T2(ii),esp));
    end
end
close;
sels = [6,16,32,48];
figure('Renderer', 'painters', 'Position', [50 50 1600 1100]);set(gcf,'DefaultLineLineWidth',2);
for ii=1:length(sels)
    subplot(1,4,ii)
    contour(T1*1e3,T2*1e3,S(:,:,sels(ii)),'LineWidth',2)
    title(['Echo #' num2str(sels(ii))])
    xlabel('T1')
    ylabel('T2')
end
sgtitle('\alpha = 120°') 
fontsize(gcf, 24, "points");
saveas(gcf,'./figs_q1/b_120.png'); close;

disp('-- Spin Echo with alpha=60 --')

alpha = 60;
num_flips = nTR*etl;
flips = pi/180*[alpha+(90-alpha/2) alpha*ones(1,num_flips-1)];
for ii=1:length(T2)
    for jj=1:length(T1)
        S(ii,jj,:) = abs(epg_cpmg(flips,[],T1(jj),T2(ii),esp));
    end
end
close;
sels = [6,16,32,48];
figure('Renderer', 'painters', 'Position', [50 50 1600 1100]);set(gcf,'DefaultLineLineWidth',2);
for ii=1:length(sels)
    subplot(1,4,ii)
    contour(T1*1e3,T2*1e3,S(:,:,sels(ii)),'LineWidth',2)
    title(['Echo #' num2str(sels(ii))])
    xlabel('T1')
    ylabel('T2')
end
sgtitle('\alpha = 60°');
fontsize(gcf, 24, "points");
saveas(gcf,'./figs_q1/b_60.png'); close;
%%