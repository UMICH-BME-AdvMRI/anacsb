%% Problem 2: SENSE
clear all;close all;clc;
location = pwd;
%% Load data
fontSize = 30;
load Data_Assignment3_Problem2.mat;
[ny,nx,nc] = size(kspaceData);
if ~exist('/figs_p2', 'dir')
       mkdir('figs_p2')
end
%% a) Fully-Sampled Image

coilCombined = sqrt(sum( coilmaps.^2, 3));
figure;
imagesc(abs(coilCombined));axis image, axis off;
colormap gray; title('Coil Sensitivity Map');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/a_Coilsens.png');
close;

allKspace = sqrt(sum(abs(kspaceData).^2,3));
figure;
imagesc(allKspace,[0 1]);axis image, axis off;
colormap gray; title('Full k-space');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/a_kspace.png');
close;

% kspace -> inverse FT -> image domain
imCombined = ifftshift(ifft2(fftshift(kspaceData)));
imCombined2 = sqrt(sum(abs(imCombined).^2,3));
figure;
imagesc(abs(imCombined2));axis image, axis off;
colormap gray; title('Magnitude image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/a_brainIm.png');
close;
orig = imCombined2;

coilIm = imCombined2.*coilCombined;
figure;
imagesc(abs(coilIm));axis image, axis off;
colormap gray; title('Coil combined Image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/a_coilComb.png');
close;
%% b) Aliased R=2 Image
R = 2;
us_ksp = zeros(ny,nx,nc);
us_ksp(R:R:ny,:,:) = kspaceData(R:R:ny,:,:);

us_kspim = sqrt(sum(abs(us_ksp).^2,3));

figure;
imagesc(abs(us_kspim),[0 1]);axis image, axis off;
colormap gray; title('Undersampled k-space');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/b_undersampkspace.png');
close;

im_us = fftshift(ifft2(ifftshift(us_ksp)));

im_us_s = sum(ifftdim(us_ksp,1:3),3);

figure;
imagesc(abs(im_us_s));axis image, axis off;
colormap gray; title('Magnitude image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/b_magimage.png');
close;

im_us_mag = sqrt(sum(abs(im_us).^2,3));
figure;
imagesc(im_us_mag);axis image, axis off;
colormap gray; title('Magnitude image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/b_magimage2.png');
close;
%% c) SENSE R=2
fov = ny/R;
im_recon = zeros(ny,nx);
sens = zeros(nc,R);
for x = 1 :nx      
    for y = 1 :fov
        % Coil image (aliased)
        I = reshape(im_us(y,x,:),nc,1);
        % Sensitivity maps per coil
        for r = 1 : R
            sens(:,r) = reshape(coilmaps(y+(r-1).*fov,x,:),nc,1);
        end
        % Pseudo inverse
        sensPinv = pinv(sens);
        rho = sensPinv*I;       % unknowns
        % SENSE image recon
        for r = 1 : R
            im_recon(y + (r-1).*fov ,x) = rho(r);
        end
    end
end

senseImage = sqrt(sum(abs(im_recon).^2,3));

figure;
imagesc(senseImage.*R);axis image, axis off;
colormap gray;title('Sense Image');colorbar();
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/c_senseR2.png');
close;

figure;
imagesc(orig - senseImage.*R);axis image, axis off;
colormap gray;title('Difference Image');colorbar();
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/c_senseDiff.png');
close;
%% d) SENCE R=4
%Aliased R=4 Image
R = 4;
us_ksp = zeros(ny,nx,nc);
us_ksp(R:R:ny,:,:) = kspaceData(R:R:ny,:,:);

us_kspim = sqrt(sum(abs(us_ksp).^2,3));

figure;
imagesc(abs(us_kspim),[0 1]);axis image, axis off;
colormap gray; title('Undersampled k-space');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/d_undersampkspace.png');
close;

im_us = fftshift(ifft2(ifftshift(us_ksp)));
im_us_s = sum(ifftdim(us_ksp,1:3),3);

figure;
imagesc(abs(im_us_s));axis image, axis off;
colormap gray; title('Magnitude Image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/d_magimage.png');
close;

im_us2 = sqrt(sum(abs(im_us).^2,3));
figure;
imagesc(abs(im_us2));axis image, axis off;
colormap gray; title('Magnitude Image');
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/d_magimage2.png');
close;

% Reconstruction
fov = ny/R;
im_recon = zeros(ny,nx);
sens = zeros(nc,R);
for x = 1 :nx      
    for y = 1 :fov
        % Coil image (aliased)
        I = reshape(im_us(y,x,:),nc,1);
        % Sensitivity maps per coil
        for r = 1 : R
            sens(:,r) = reshape(coilmaps(y+(r-1).*fov,x,:),nc,1);
        end
        % Pseudo inverse
        sensPinv = pinv(sens);
        rho = sensPinv*I;
        % SENSE image recon            
        for r = 1 : R
            im_recon(y + (r-1).*fov ,x) = rho(r);
        end
    end
end

senseImage = sqrt(sum(abs(im_recon).^2,3));
senseImage = senseImage .* R;  % to match intensities with original

figure;
imagesc(senseImage);axis image, axis off;
colormap gray;title('Sense Image');colorbar();
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/d_senseR4.png');
close;

figure;
imagesc(orig - senseImage);axis image, axis off;
colormap gray;title('Difference Image');colorbar();
fontsize(gcf, fontSize, "points");
saveas(gcf,'./figs_p2/d_senseDiff.png');
close;
