%% Problem 1: Partial Fourier Imaging
clear all;close all;clc;
location = pwd;
%% Load data
fontSize = 30;
load Data_Assignment3_Problem1.mat;
if ~exist('/figs_p1', 'dir')
       mkdir('figs_p1')
end
%% a) Zero-Filled Recon
[ny,nx] = size(kspaceData_SingleCoil);
% Plot
figure;
imagesc(abs(kspaceData_SingleCoil),[0 1]);axis image, axis off;
colormap gray; title('Full k-space');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p1/a_fullKspace.png');
close;

figure;
imagesc(abs(ifftdim(kspaceData_SingleCoil,1:2)));axis image, axis off;
colormap gray; title('Image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p1/a_fullImage.png');
close;

% zero-filling
factor = 5/8;
ky = factor*ny;
kspace_zf = zeros(ny,nx);
kspace_zf(1:ky,:) = kspaceData_SingleCoil(1:ky,:);
im_zf = ifftdim(kspace_zf,1:2);
figure;
imagesc(abs(kspace_zf),[0 1]);axis image, axis off;
colormap gray; title('Partial k-space');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p1/a_partialKspace.png');
close;

figure;
subplot(121);
imagesc(abs(im_zf));colormap gray; title('Magnitude');axis image, axis off;
fontsize(gcf, fontSize, "points");
subplot(122);
imagesc(angle(im_zf),[-pi pi]);axis image, axis off;
colormap gray; title('Phase');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p1/a_partialKspaceMagPhase.png');
close;

% full-kspace
im_recon = ifftdim(kspaceData_SingleCoil,1:2);
figure,
subplot(121);
imagesc(abs(im_recon));axis image, axis off;
colormap gray; title('Magnitude');axis image, axis off;
fontsize(gcf, fontSize, "points");
subplot(122);
%imagesc(angle(im_recon).*(abs(im_recon)>1e-3),[-pi pi]);axis image, axis off;
imagesc(angle(im_recon),[-pi pi]);axis image, axis off;
colormap gray; title('Phase');
fontsize(gcf, fontSize, "points"); 
close;

%difference
im_diff = im_recon - im_zf;
lm = abs(im_diff(:));
figure,
subplot(121);
imagesc(abs(im_diff),[0+1e-5 max(lm)-0.0005]);axis image, axis off;
colormap gray; title('Mag. diff');
fontsize(gcf, fontSize, "points");
subplot(122);
imagesc(angle(im_diff),[-pi pi]);axis image, axis off;
colormap gray; title('Phase diff');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p1/a_Diff.png');
close;
%% b) POCS

% Step 1: Phase estimation
center = ny/2;
numLines = ky - center;
numIni = center - numLines;
numEnd = center + numLines;
indx = numIni:numEnd;

win = hann(numel(indx), 'symmetric' );
winAll = zeros(ny,1);
winAll(indx) = win;

filter  = repmat(winAll,[1 200]);
ks_phase = kspace_zf.*filter;

figure;
imagesc(abs(ks_phase),[0 1]); colormap gray;
title('kspace multiplied by Hann filter');
fontsize(gcf, fontSize, "points"); close;

im_pksp = ifftdim(kspaceData_SingleCoil,1:2);
imMag = ifftshift(ifft2(fftshift(ks_phase)));
imPhase = angle(imMag);
kspCorrected = kspace_zf;

iter = 1:10;

for ii = iter
    if ii == 1
        im_pksp = ifftdim(kspace_zf,1:2);
    else
        im_pksp = ifftdim(ks_phase,1:2);
    end

    imNew = abs(im_pksp).*exp(-i*imPhase);
    kspNew = fftshift(fft2(ifftshift(imNew)));
    
    kspCorrected(numEnd+1:end,:) = kspNew(numEnd+1:end,:);
    ks_phase = kspCorrected;
    if ii == max(iter)
        im_pksp = ifftdim(ks_phase,1:2);
    end
end

% Step 1:
figure('Renderer', 'painters', 'Position', [50 50 1000 500]);set(gcf,'DefaultLineLineWidth',2);subplot(121);
imagesc(abs(imMag)); colormap gray;axis image, axis off;
title('Magnitude');
fontsize(gcf, fontSize, "points");
subplot(122);
imagesc(abs(imPhase)); colormap gray;axis image, axis off;
title('Phase');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p1/b_step1.png');
close;

% Step 2: Iterative Reconstruction
figure('Renderer', 'painters', 'Position', [50 50 1000 500]);set(gcf,'DefaultLineLineWidth',2);subplot(121);
subplot(121);
imagesc(abs(im_pksp));axis image, axis off;
colormap gray; title('POCS - Magnitude');axis image, axis off;
fontsize(gcf, fontSize, "points");
subplot(122);
imagesc(angle(im_pksp),[-pi pi]);axis image, axis off;
colormap gray; title('POCS - Phase');
fontsize(gcf, fontSize, "points"); 
saveas(gcf,'./figs_p1/b_POCS.png');
close;

%difference
im_diff = im_recon - im_pksp;
lm = abs(im_diff(:));
figure,
subplot(121);
imagesc(abs(im_diff),[0+1e-5 max(lm)-0.0005]);axis image, axis off;
colormap gray; title('Mag. diff');
fontsize(gcf, fontSize, "points");
subplot(122);
imagesc(angle(im_diff),[-pi pi]);axis image, axis off;
colormap gray; title('Phase diff');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p1/b_diffPOCS.png');
close;
