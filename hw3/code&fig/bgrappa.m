%% GRAPPA
clear all;close all;clc;
location = pwd;
%% Load data
fontSize = 30;
load Data_Assignment3_Problem2.mat;
[ny,nx,nc] = size(kspaceData);
if ~exist('/figs_bns', 'dir')
       mkdir('figs_bns')
end
originalImage = ifftshift(ifft2(fftshift(kspaceData)));
originalImage = sqrt(sum(abs(originalImage).^2,3));
%% full and undersampled kspace
% original data, full kspace
brainData = ifftshift(ifft2(fftshift(kspaceData)));
im_us_mag = sqrt(sum(abs(brainData).^2,3));

% full kspace
allKspace = sqrt(sum(abs(kspaceData).^2,3));
figure;
imagesc(allKspace,[0 1]);axis image, axis off;
colormap gray; title('Full k-space');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_bns/fullKspace.png');
close;

% undersample ks-ace
% 3(readout) x 2(phase encoding) -> undersample R=2
R = 2;
us_ksp = zeros(ny,nx,nc);
us_ksp(1:R:ny,:,:) = kspaceData(1:R:ny,:,:);

us_kspim = sqrt(sum(abs(us_ksp).^2,3));

figure;
imagesc(abs(us_kspim),[0 1]);axis image, axis off;
colormap gray; title('Undersampled k-space');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_bns/usKspace.png');
close;

usImage = ifftshift(ifft2(fftshift(us_ksp)));
us_im = sqrt(sum(abs(usImage).^2,3));

figure;
imagesc(us_im);axis image, axis off;
colormap gray; title('Undersampled Image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_bns/usImage.png');
close;

%% ACS
acsLines = 24;
kernelRO = 3;
kernelPE = 2;
kernelSize = kernelPE .* kernelRO .* nc;

acsIni = ny/2 - acsLines/2;
ACS = zeros(ny,nx,nc);
ACS(acsIni:acsIni+acsLines-1,:,:) = kspaceData(acsIni:acsIni+acsLines-1,:,:);
ACSKspace = sqrt(sum(abs(ACS).^2,3));
figure;
imagesc(ACSKspace,[0 1]);axis image, axis off;
colormap gray; title('ACS');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_bns/ACS.png');
close;
%% Source and targets in ACS ---> find weights
finder = any(us_ksp, [2, 3]);
% Source
kNo = 1; % kernel/patch number
for ky = acsIni : acsIni+acsLines-1
    for kx = 1:(ny-R)  

        for kk = 1 :nc
            tmp(:,kk) = [ACS(ky,kx:kx+R,kk).'; ACS(ky+R,kx:kx+R,kk).'];
        end

        S(kNo,:) = reshape(tmp,[kernelSize,1]); % ch1: each patch 3x3 is reshaped into vector and put into matrix one line after another
       
        kNo = kNo + 1; % to move through all patches
        clear tmp
    end
end
%Targets
kNo = 1; % kernel/patch number
for ky = acsIni : acsIni+acsLines-1
    for kx = 1:(ny-R)  

        tmp = ACS(ky+1,kx+1,:);

        T(kNo,:) = reshape(tmp,[nc,1]); % ch1: each patch 3x3 is reshaped into vector and put into matrix one line after another
       
        kNo = kNo + 1; % to move through all patches
        clear tmp
    end
end

W = pinv(S)*T;
%% Find missing lines
kNo = 1; % kernel/patch number
imrec = us_ksp;
for ky = 1:R:(ny-R)   %N-2, N-3 for R=4
    for kx = 1:(nx-R)
        for kk = 1 : nc
            tmp(:,kk) = [us_ksp(ky,kx:kx+R,kk).'; us_ksp(ky+R,kx:kx+R,kk).'];
        end

        S_new(kNo,:) = reshape(tmp,[1,kernelSize]);
       
        targets = S_new(kNo,:)*W;

        imrec(ky+1,kx+1,:)=targets;
        kNo = kNo + 1;
        clear tmp
    end
end

%T_new = S_new*W;
%% Filling missing lines
% Ny = length(1:R:(ny-R));
% Nx = length(1:(nx-R));
% TTargets = reshape(T_new,[Nx,Ny,nc]);
% us_kspnew = us_ksp;
% for kk = 1 :nc
%     tmp = TTargets(:,:,kk);
%     us_kspnew(R:R:ny-1,1:nx-R,kk) = tmp.'; 
% end
us_kspn = imrec;
us_kspn(acsIni:acsIni+acsLines-1,:,:) = kspaceData(acsIni:acsIni+acsLines-1,:,:);

%us_kspimn = sqrt(sum(abs(us_kspn).^2,3));

for kk = 1 :nc
    tmp = us_kspn(R:end-1,R:end-1,kk);
    imGrappa(:,:,kk) = fftshift(ifft2(ifftshift(tmp)));
end

grappa_recon = sqrt(sum(abs(imGrappa).^2,3));
figure;
imagesc(grappa_recon);axis image, axis off;
colormap gray;title('GRAPPA Image');colorbar();
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_bns/grappaImage.png');
close;

originalImage = originalImage(R:end-1,R:end-1);
figure;
imagesc(originalImage - grappa_recon);axis image, axis off;
colormap gray;title('Difference Image');colorbar();
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_bns/grappaDif.png');
close;