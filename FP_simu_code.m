close all;clear;clc;

%% The high resolution complex object
% Add necessary folders into the current working directory
addpath(genpath(pwd));

% Load data file
data_name = 't000';
data_dir = ['Data\' 'DIC-C2DH-HeLa\' data_name '.tif']; 

objectAmplitude = read(Tiff(data_dir));
% phase = double(imread('westconcordorthophoto.png'));
% phase = pi*imresize(phase,[256 256])./max(max(phase));
object = mat2gray(objectAmplitude,[50 200]);
figure;
imshow(object,[]);
title('HR image');

%% setup the parameters of simulate/experiment system
wlength = 6.3000e-07; % illu wavelength, in m
z = 0; % defocus distance, in m
aberration = 0; % pre-calibrate aberration
xint = 0;yint = 0; % offset of initial LED to the patch center, in mm 
theta = 0; % rotation angle of LED array to the camera sensor frame, in degree
xstart = 18; ystart = 20; % absolute coordinate of initial LED
arraysize = 15; % side length of lit LED array
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
H      = 90.88; % distance between LEDs and sample, in mm
LEDp   = 4;     % distance between adjacent LEDs, in mm
nglass = 1.52;  % refraction index of glass substrate
t      = 1;     % glass thickness, in mm
[kx, ky, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize^2);

NA          = 0.1;      % objective NA
spsize      = 1.845e-6; % pixel size of low-res image on sample plane, in m
upsmp_ratio = 4;        % upsampling ratio
psize       = spsize/upsmp_ratio; % pixel size of high-res image on sample plane, in m

%% k-space parameterization
[m, n] = size(object);
pratio = round(spsize/psize); % upsampling ratio
m1 = m/pratio; n1 = n/pratio;
k0 = 2*pi/wlength;
kx = k0*kx; ky = k0*ky;
NAfilx = NA*(1/wlength)*n*psize; NAfily = NA*(1/wlength)*m*psize; % m1*spsize = m*psize
kmax = pi/psize; % the max wave vector of the OTF
dkx = 2*pi/(psize*n); dky = 2*pi/(psize*m);
kx2 = -kmax:kmax/((n-1)/2):kmax; ky2 = -kmax:kmax/((m-1)/2):kmax; % odd N
[kxm, kym] = meshgrid(kx2,ky2); kzm = sqrt(k0^2-kxm.^2-kym.^2);

%% prior knowledge of aberration
H2 = exp(1j.*z.*real(kzm)).*exp(-abs(z).*abs(imag(kzm))); % define the defocus aberration if it is known or you want to test it
astigx = 0; astigy = 0; % define the astigmatism aberration if it is known or you want to test it
[M1, N1] = meshgrid(1:m1,1:n1);
zn = astigx*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),2,2)+...
     astigy*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),-2,2);
zn = imresize(zn,[m1,n1]);
if aberration ~= 0
    fmaskpro = aberration; % pre-calibrated aberrations
else
    fmaskpro = 1.*double(((N1-(m1+1)/2)/NAfily).^2+((M1-(n1+1)/2)/NAfilx).^2<=1)... % low-pass filter
    .*H2(round((m+1)/2-(m1-1)/2):round((m+1)/2+(m1-1)/2),round((n+1)/2-(n1-1)/2):round((n+1)/2+(n1-1)/2))... % defocus aberration
    .*exp(pi*1j.*zn); % astigmatism aberration
    % In this example, we can test the effect of the defocus aberration (z) and astigmatism aberrations (astigx, astigy) 
    % If the aberration is unknown, one can test different z, astigx, astigy for the best result
    % Gradient descent can also be used to update z, astigx, astigy (did not implement here)
    % Higher-order Zernike modes can be tested in a similar manner
end

%% forward image processing
objectFT = fftshift(fft2(object));
imlow_HDR = zeros(m1,n1,arraysize^2);
for i3 = 1:arraysize^2        
        kxc=round((n+1)/2-kx(1,i3)/dkx);
        kyc=round((m+1)/2-ky(1,i3)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        O_j=objectFT(kyl:kyh,kxl:kxh);
        imgseqlowFT=O_j.*fmaskpro;
        imlow_HDR(:,:,i3)=abs(ifft2(ifftshift(imgseqlowFT))); 
end


figure;
set(gcf,'outerposition',get(0,'ScreenSize'))
imshow(imlow_HDR(:,:,1),[]);
title(['LR image ' num2str(1)]);

%% save the simulate result 
save('Data\HeLa.mat', 'aberration', 'imlow_HDR', 'theta', 'wlength', 'xint', 'yint', 'z');