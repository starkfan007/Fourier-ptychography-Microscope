function seq = sparse_recover(imseqlow,kx,ky,NA,wlength,spsize,psize,z, opts)
% FP algorithm to recover high-resolution image from low-resolution measured images
% Input:
%       imseqlow: low-res measurements, [m1 x n1 x numim] matrix
%       kx,ky: normalized wavevector of each LED, which should times k0 later
%       NA: objective NA
%       wlength: central peak wavelength of LED, in m
%       spsize: pixel size of low-res image on sample plane, in m
%       psize: pixel size of high-res image on sample plane, in m
%       z: known defocus distance, in m
%       opts: other optimization parameters
% Output:
%       LED_idx: LED index that are used to recover HR Image
%       imseqlow_FT: low-res Image freq spectrum are used truly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(opts,'aberration')
    opts.aberration = 0;
end
aberration = opts.aberration;

%% k-space parameterization
[m1, n1, numim] = size(imseqlow);
pratio = round(spsize/psize); % upsampling ratio
m = pratio*m1; n = pratio*n1;
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

%% initialization
him = imresize(imseqlow(:,:,1),[m,n]); 
him = mat2gray(him);
himFT = fftshift(fft2(him));
imgseqlowlow = zeros(m1, n1, numim);

%% find the main spectrum region     
for i3 = 1:numim         
    kxc=round((n+1)/2-kx(1,i3)/dkx);
    kyc=round((m+1)/2-ky(1,i3)/dky);
    kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
    kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
    O_j=himFT(kyl:kyh,kxl:kxh);
    lowFT=O_j.*fmaskpro;
    im_lowFT=ifft2(ifftshift(lowFT));
    imgseqlowlow(:,:,i3)=abs(im_lowFT);
end
seq = calc_entropy(imgseqlowlow);

end

