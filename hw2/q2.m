%% Problem 2: Single and MSE sequences
clear all;close all;clc;
addpath('functions')
load('brain_maps.mat');
%% a) Single-echo spin echo
% Min T1 and T2 values
fprintf('Minimun T1 value %d ms\n',min(nonzeros(T1map(:))));
fprintf('Maximum T1 value %d ms\n',max(nonzeros(T1map(:))));
fprintf('Av. T1 value %d ms\n',round(mean(nonzeros(T1map(:)))));

fprintf('Minimun T2 value %d ms\n',min(nonzeros(T2map(:))));
fprintf('Maximum T2 value %d ms\n',max(nonzeros(T2map(:))));
fprintf('Av. T2 value %d ms\n',round(mean(nonzeros(T2map(:)))));

% min. T1 is 10ms
% min. T2 is 2 ms
%   (1 - e^-(TR/T1))*e^(-TE/T2)

% T1w -> short TR and short TE -> TE < 79 ms (av. value) and TR < 777 ms
% T2w -> long TR and long TE -> TR > 777 ms and TE > 79 ms
% PD -> long TR, short TE -> TR > 777 ms and TE < 79 ms
%% T1w
TE = 15; TR = 500;        % ms
dT = 1;		                % 1ms delta-time.
df = 0;

for xx = 1 :  size(T1map,1)
    for yy = 1 :  size(T1map,2)
        [Msig,~] = sesignal(T1map(xx,yy),T2map(xx,yy),TE,TR,0);
        T1w(xx,yy) = Msig;
    end
end
imagesc(abs(T1w));axis off;title(sprintf('T1w, TE = %d ms',TE));colormap gray;
if ~exist('/figs_q2', 'dir')
       mkdir('figs_q2')
end
saveas(gcf,'./figs_q2/a_T1w.png');
fontsize(gcf, 30, "points"); 
close;
%% PD
TE = 15; TR = 4000;        % ms
dT = 1;		                % 1ms delta-time.
df = 0;

for xx = 1 :  size(T1map,1)
    for yy = 1 :  size(T1map,2)
        [Msig,~] = sesignal(T1map(xx,yy),T2map(xx,yy),TE,TR,0);
        PD(xx,yy) = Msig;
    end
end
imagesc(abs(PD));axis off;title(sprintf('PD, TE = %d ms',TE));colormap gray;
if ~exist('/figs_q2', 'dir')
       mkdir('figs_q2')
end
saveas(gcf,'./figs_q2/a_PD.png');
fontsize(gcf, 30, "points"); 
close;
%% T2w
TE = 100; TR = 6000;        % ms
dT = 1;		                % 1ms delta-time.
df = 0;

for xx = 1 :  size(T1map,1)
    for yy = 1 :  size(T1map,2)
        [Msig,~] = sesignal(T1map(xx,yy),T2map(xx,yy),TE,TR,0);
        T2w(xx,yy) = Msig;
    end
end
imagesc(abs(T2w));axis off;title(sprintf('T2w, TE = %d ms',TE));colormap gray;
if ~exist('/figs_q2', 'dir')
       mkdir('figs_q2')
end
saveas(gcf,'./figs_q2/a_T2w.png');
fontsize(gcf, 30, "points"); 
close;

%% b) Fast spin echo
%% i)
clear all; close all; clc;
nTR = 5;
esp = 5/1000;      % sec
etl = 32;         % echo train length
TR = 3*1e3;       % ms
TE = 100;         % ms
N = 2*etl;
Qall = zeros(4*etl,etl,nTR);    % store f, F* states
Zall = zeros(2*etl,etl,nTR);    % store Z states
Q = epg_m0(N); 
T1 = [1000 1000 2000 2000];              % ms
T2 = [50 100 50 100];                 % ms
alpha = 180;
num_flips = 32;
flips = pi/180*[alpha+(90-alpha/2) alpha*ones(1,num_flips-1)];

nTR = 5; 
kk = 0;

for jj = 1: nTR
    S(1 + (etl*kk) : etl*(kk + 1)) = abs(epg_cpmg(flips,etl,T1(1)/1000,T2(1)/1000,esp));
    kk = kk +1;
end

for nn = 1 :  length(T1)
    kk = 0;
    for tt = 1 : nTR
        [Msig,Mss] = fsesignal2(T1(nn),T2(nn),TE,TR,0,etl);
        signal(nn,1 + (etl*kk) : etl*(kk + 1)) = Msig; 
        kk = kk +1;
    end
end

timet = (1:32)*esp+12;
figure('Renderer', 'painters', 'Position', [50 50 1200 1000]);set(gcf,'DefaultLineLineWidth',2);
legend1 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(1),T2(1));
legend2 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(2),T2(2));
legend3 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(3),T2(3));
legend4 = sprintf('T_{1} = %d ms, T_{2} = %d ms',T1(4),T2(4));
plot(timet,abs(signal(:,129:end))); title('Transverse Magnetization vs echo time');
legend(legend1,legend2,legend3,legend4)
lplot('Echo Time (ms)','Signal','Signal vs. Echo Time');
fontsize(gcf, 28, "points");
saveas(gcf,'./figs_q2/b_i.png');
close;

%% ii) FSE to create an image
clear all;close all;clc;
location = '/Users/anacecisb/Documents/MATLAB/BME599/hw2';
cd(location)
load('brain_maps.mat');
etl = 32;         % echo train length
TR = 3*1e3;           % sec
TE =  5;              % ms
[Nx , Ny] = size(T1map);
nTR = Ny/etl;      %number of TR

fprintf('Total acquisition time: %.2f s\n', nTR*TR/(1e3));

alpha = pi;
flipangle = [(alpha + (pi/2 - alpha/2))  ones(1,etl-1)*alpha];
ksp = zeros(Nx,Ny);
signal = zeros(size(Nx,1),size(Ny,2),etl*nTR);
for ky = 1 : Ny
    fprintf('Line %d\n',ky);
    for kx = 1 : Nx
            kk = 0;
        for jj = 1: nTR
            [Msig,Mss] = fsesignal(T1map(kx,ky),T2map(kx,ky),TE,TR,0,etl,M0map(kx,ky));
            signal(kx,ky,1 + (etl*kk) : etl*(kk + 1)) = Msig; 
            kk = kk +1;
        end
    end
end

%% kspace
esp = 5*1e-3;
TEeff = 80;         %ms
nespEff = (TEeff*1e-3)/esp;

blocks = Ny/etl;

for bb = 1 :  blocks
    inds(:,bb) = [Ny+1-bb-nTR:-nTR:1  Ny+1-bb];
end

% --------testing 32 blocks
inds2 = inds(:);
ksp = zeros(Nx,Ny);
for ech = 1 :  Ny            % all kspace lines y-axis
    tmp = signal(:,:,ech);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    ksp(inds2(ech),:) = tmp2(inds2(ech),:);
    %imagesc(abs(ksp));pause(0.1)
end

figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
imagesc(abs(ksp));colormap gray; axis off;title('kspace');
fontsize(gcf, 30, "points");
saveas(gcf,'./figs_q2/b_ii_kspace.png');close;

% Recovering the image
img = (ifft2(ifftshift(ksp)));
figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
imagesc(abs(img));colormap gray; axis off;title('Image recon');
fontsize(gcf, 30, "points");
saveas(gcf,'./figs_q2/b_ii_img.png');close;

fprintf('FSE total scan time: %.2f sec \n',nTR*TR);
fprintf('Single-echo spin echo total scan time: %.2f sec \n',Ny*TR);
%% iii)
esp = 5*1e-3;
TEeff = 40;         %ms
nespEff = (TEeff*1e-3)/esp;
Neff = 128;
nnTEini = [Neff-nespEff+1:Neff-nespEff+etl Neff-nespEff:-1:Neff-nespEff+1-etl Neff-nespEff+etl+1:Neff-nespEff+2*etl...
            Neff-nespEff-etl:-1:Neff-nespEff+1-2*etl Neff-nespEff+2*etl+1:Neff-nespEff+3*etl...
            Neff-nespEff-2*etl:-1:Neff-nespEff+1-3*etl...
            Neff-nespEff+3*etl+1:Neff-nespEff+4*etl Neff-nespEff-3*etl:-1:1....
            Neff-nespEff+4*etl+1:Neff*2];

ksp40 = zeros(Nx,Ny);
for ech = 1 :  Ny            % all kspace lines y-axis
    tmp = signal(:,:,ech);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    ksp40(nnTEini(ech),:) = tmp2(nnTEini(ech),:);
    %imagesc(abs(ksp40));pause(0.1)
end

figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
imagesc(abs(ksp40));colormap gray; axis off;title('kspace, TE_{eff} = 40 ms');
fontsize(gcf, 30, "points");
saveas(gcf,'./figs_q2/b_ii_kspace40.png');close;

% Recovering the image
img40 = (ifft2(ifftshift(ksp40)));
figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
imagesc(abs(img40));colormap gray; axis off;title('Image recon, TE_{eff} = 40 ms');
fontsize(gcf, 30, "points");
saveas(gcf,'./figs_q2/b_ii_img40.png');close;

%----------- TEeff = 120
TEeff = 120;         %ms
nespEff = (TEeff*1e-3)/esp;
Neff = 128;
clear inds
for bb = 1 :  blocks/2
    inds(:,bb) = [220-bb+1:-4:1 256-bb+1:-4:221]%[Ny+1-bb:-nTR:1];
end

nnTEini = inds(:);
for ky = 1 :  Ny            % all kspace lines y-axis
    tmp = signal(:,:,ky);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    ksp120(nnTEini(ky),:) = tmp2(nnTEini(ky),:);
end

figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
imagesc(abs(ksp120));colormap gray; axis off;title('kspace, TE_{eff} = 120 ms');
fontsize(gcf, 30, "points");
saveas(gcf,'./figs_q2/b_ii_kspace120.png');close;

% Recovering the image
img120 = (ifft2(ifftshift(ksp120)));
figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
imagesc(abs(img120));colormap gray; axis off;title('Image recon, TE_{eff} = 120 ms');
fontsize(gcf, 30, "points");
saveas(gcf,'./figs_q2/b_ii_img120.png');close;

%% iv)
% simulation at different ETL
clear all;close all;clc;
location = '/Users/anacecisb/Documents/MATLAB/BME599/hw2';
cd(location)
load('brain_maps.mat');
TR = 3*1e3;           % sec
TE =  5;              % ms
[Nx , Ny] = size(T1map);
step = Ny/32;
etl = 16;         % echo train length
nTR = Ny/etl;      %number of TR
signaletl = zeros(size(Nx,1),size(Ny,2),etl*nTR);
disp('---- ETL = 16 ----');
for ky = 1 : Ny
    fprintf('Line %d\n',ky);
    for kx = 1 : Nx
            kk = 0;
        for jj = 1: nTR
            [Msig,Mss] = fsesignal(T1map(kx,ky),T2map(kx,ky),TE,TR,0,etl,M0map(kx,ky));
            signaletl(kx,ky,1 + (etl*kk) : etl*(kk + 1)) = Msig; 
            kk = kk +1;
        end
    end
end

esp = 5*1e-3;
TEeff = 80;         %ms
nespEff = (TEeff*1e-3)/esp;

blocks = Ny/etl;
for bb = 1 :  blocks/2
    inds(:,bb) = [Ny+1-bb-step:-step:1  Ny+1-bb];
end

inds2 = inds(:);
kspetl = zeros(Nx,Ny);
for ech = 1 :  Ny            % all kspace lines y-axis
    tmp = signaletl(:,:,ech);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    kspetl(inds2(ech),:) = tmp2(inds2(ech),:);
    %imagesc(abs(kspetl));pause(0.1)
end

    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(kspetl));colormap gray; axis off;title(sprintf('kspace, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_kspetl_16.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

% Recovering the image
    imgetl= ifft2(ifftshift(kspetl(:,:)));
    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(imgetl));colormap gray; axis off;title(sprintf('Image recon, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_imgetl_16_v2.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;


disp('---- ETL = 32 ----');
etl = 32;         % echo train length
nTR = Ny/etl;      %number of TR

signaletl = zeros(size(Nx,1),size(Ny,2),etl*nTR);
for ky = 1 : Ny
    fprintf('Line %d\n',ky);
    for kx = 1 : Nx
            kk = 0;
        for jj = 1: nTR
            [Msig,Mss] = fsesignal(T1map(kx,ky),T2map(kx,ky),TE,TR,0,etl,M0map(kx,ky));
            signaletl(kx,ky,1 + (etl*kk) : etl*(kk + 1)) = Msig; 
            kk = kk +1;
        end
    end
end

esp = 5*1e-3;
TEeff = 80;         %ms
nespEff = (TEeff*1e-3)/esp;
clear inds
blocks = Ny/etl;
for bb = 1 :  blocks
    inds(:,bb) = [Ny+1-bb-nTR:-nTR:1  Ny+1-bb];
end

inds2 = inds(:);
kspetl = zeros(Nx,Ny);
for ech = 1 :  Ny            % all kspace lines y-axis
    tmp = signaletl(:,:,ech);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    kspetl(inds2(ech),:) = tmp2(inds2(ech),:);
   % imagesc(abs(kspetl));pause(0.1)
end

    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(kspetl));colormap gray; axis off;title(sprintf('kspace, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_kspetl_32.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

% Recovering the image
    imgetl= ifft2(ifftshift(kspetl(:,:)));
    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(imgetl));colormap gray; axis off;title(sprintf('Image recon, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_imgetl_32.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;


    disp('---- ETL = 64 ----');
etl = 64;         % echo train length
nTR = Ny/etl;      %number of TR
step = Ny/32;
signaletl = zeros(size(Nx,1),size(Ny,2),etl*nTR);
for ky = 1 : Ny
    fprintf('Line %d\n',ky);
    for kx = 1 : Nx
            kk = 0;
        for jj = 1: nTR
            [Msig,Mss] = fsesignal(T1map(kx,ky),T2map(kx,ky),TE,TR,0,etl,M0map(kx,ky));
            signaletl(kx,ky,1 + (etl*kk) : etl*(kk + 1)) = Msig; 
            kk = kk +1;
        end
    end
end

esp = 5*1e-3;
TEeff = 80;         %ms
nespEff = (TEeff*1e-3)/esp;
blocks = Ny/etl;
clear inds
for bb = 1 :  blocks
    inds(:,bb) = [Ny+1-bb-(nTR):-nTR:1  Ny+1-bb];
end

inds2 = inds(:);
kspetl = zeros(Nx,Ny);
for ech = 1 :  Ny            % all kspace lines y-axis
    tmp = signaletl(:,:,ech);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    kspetl(inds2(ech),:) = tmp2(inds2(ech),:);
    %imagesc(abs(kspetl));pause(0.1)
end

    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(kspetl));colormap gray; axis off;title(sprintf('kspace, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_kspetl_64.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

% Recovering the image
    imgetl= ifft2(ifftshift(kspetl(:,:)));
    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(imgetl));colormap gray; axis off;title(sprintf('Image recon, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_imgetl_64.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

    disp('---- ETL = 128 ----');
etl = 128;         % echo train length
nTR = Ny/etl;      %number of TR

signaletl = zeros(size(Nx,1),size(Ny,2),etl*nTR);
for ky = 1 : Ny
    fprintf('Line %d\n',ky);
    for kx = 1 : Nx
            kk = 0;
        for jj = 1: nTR
            [Msig,Mss] = fsesignal(T1map(kx,ky),T2map(kx,ky),TE,TR,0,etl,M0map(kx,ky));
            signaletl(kx,ky,1 + (etl*kk) : etl*(kk + 1)) = Msig; 
            kk = kk +1;
        end
    end
end

esp = 5*1e-3;
TEeff = 80;         %ms
nespEff = (TEeff*1e-3)/esp;
clear inds
blocks = Ny/etl;
for bb = 1 :  blocks
    inds(:,bb) = [Ny+1-bb-(nTR):-nTR:1  Ny+1-bb];
end

inds2 = inds(:);
kspetl = zeros(Nx,Ny);
for ech = 1 :  Ny            % all kspace lines y-axis
    tmp = signaletl(:,:,ech);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    kspetl(inds2(ech),:) = tmp2(inds2(ech),:);
   % imagesc(abs(kspetl));pause(0.1)
end

    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(kspetl));colormap gray; axis off;title(sprintf('kspace, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_kspetl_128.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

% Recovering the image
    imgetl= ifft2(ifftshift(kspetl(:,:)));
    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(imgetl));colormap gray; axis off;title(sprintf('Image recon, TE_{eff} = 80 ms, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_imgetl_128.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;
% --------- Single spin echo image with TE = 80 ms

tmp = signaletl(:,:,nespEff);
tmp(isnan(tmp))=0;   
ksp = fftshift(fft2(tmp));
 figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(ksp));colormap gray; axis off;title(sprintf('kspace, TE_{eff} = 80 ms, ETL = %d, Single spin echo',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_kspetl_SE80ms.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

   img= ifft2(ifftshift(ksp));
    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    imagesc(abs(img));colormap gray; axis off;title(sprintf('Image recon, TE_{eff} = 80 ms, Single spin echo',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_ii_imgetl_SE80ms.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;


%% c) bonus problem
clear all;close all;clc;
disp('---- ETL = 128 ----');
location = '/Users/anacecisb/Documents/MATLAB/BME599/hw2';
cd(location)
load('brain_maps.mat');
TR = 3*1e3;           % ms
TE =  5;              % ms
[Nx , Ny] = size(T1map);
etl = 128;         % echo train length
nTR = Ny/etl;      %number of TR
step = Ny/32;

signaletl = zeros(size(Nx,1),size(Ny,2),etl*nTR);
for ky = 1 : Ny
    fprintf('Line %d\n',ky);
    for kx = 1 : Nx
            kk = 0;
        for jj = 1: nTR
            [Msig,Mss] = fsesignal(T1map(kx,ky),T2map(kx,ky),TE,TR,0,etl,M0map(kx,ky));
            signaletl(kx,ky,1 + (etl*kk) : etl*(kk + 1)) = Msig; 
            kk = kk +1;
        end
    end
end
esp = 5*1e-3;
TEeff = 80;         %ms
nespEff = (TEeff*1e-3)/esp;
Neff = 128;

nnTEini = [Neff-nespEff+1:Neff-nespEff+etl 112:-1:1 241:256];

hzeros = 112;               % (256*7)/16 zeros
inds2 = nnTEini(:);
kspetl = zeros(Nx,Ny);
for ech = 1 :  Ny            % all kspace lines y-axis
    tmp = signaletl(:,:,ech);
    tmp(isnan(tmp))=0;   
    tmp2 = fftshift(fft2(tmp));
    kspetl(inds2(ech),:) = tmp2(inds2(ech),:);
   %imagesc(abs(kspetl));pause(0.1)
end

kspdata = kspetl;
kspdata(1+Nx-hzeros:end,:) = 0;
img_zf = abs(ifftdim(kspdata,1:2));
imz = fftshift(ifftn(fftshift(kspdata)));

data_center = kspdata;
data_center(1:hzeros,:) = 0;
im_ph = fftshift(ifftn(fftshift(data_center)));
im_pc = imz.*exp(-1i*angle(im_ph));

datapc = fftshift(fftn(fftshift(im_pc)));
datapc(1+Nx-hzeros:end,:) = 0;

datapc(1+Nx-hzeros:end,:) = rot90(datapc(1:hzeros,:),2);

impc = (ifft2(ifftshift(datapc)));

    figure('Renderer', 'painters', 'Position', [50 50 1100 900]);set(gcf,'DefaultLineLineWidth',2);
    subplot(121)
    imagesc(log10(abs(kspdata)));colormap gray; axis off;title(sprintf('Partial kspace, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    subplot(122)
    imagesc(log10(abs(datapc)));colormap gray; axis off;title(sprintf('Full kspace, ETL = %d',etl));
    fontsize(gcf, 30, "points");
    fname = ['b_c_ksp.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

% Recovering the image, conjugate
    figure
    imagesc(abs(impc));colormap gray; axis off;title('Image recon, conjugate synthesis');
    fontsize(gcf, 30, "points");
    fname = ['b_c_img.png'];
    saveas(gcf,strcat('./figs_q2/',fname));close;

    % Pocs
    imReco = pocs( kspdata, 30, true );

close;
figure,
imagesc(abs(ifftshift(imReco)));
axis off; colormap gray;
title(' POCS reco');fontsize(gcf, 30, "points");
fname = ['b_c_pocs.png'];
saveas(gcf,strcat('./figs_q2/',fname));close;
% Zero
figure
    imagesc(abs(ifftshift(img_zf)));colormap gray; axis off;title('Image recon, zero filling');
    fontsize(gcf, 30, "points");
 fname = ['b_c_zero.png'];
saveas(gcf,strcat('./figs_q2/',fname));close;