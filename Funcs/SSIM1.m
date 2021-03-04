function [mssim, ssim_map,siga_sq,sigb_sq] = SSIM1(ima, imb)  
% ========================================================================    
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image  
% quality assessment: From error visibility to structural similarity,"  
% IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612,  
% Apr. 2004.  
% Add window: w=fspecial('gaussian', 11, 1.5);  
%                 (2*ua*ub+C1)*(2*sigmaa*sigmab+C2)  
%   SSIM(A,B)=！！！！！！！！！！！！！！！！！！！！！！！！  
%              (ua*ua+ub*ub+C1)(sigmaa*sigmaa+sigmab*sigmab+C2)  
%     C1=(K1*L);  
%     C2=(K2*L);   K1=0.01,K2=0.03  
%     L is gray scale, L=255  
%-------------------------------------------------------------------  
%     ima - img A  
%     imb - img B  
%  
% ssim_map - SSIM(A,B|w)  
%    mssim - average SSIM(A,B|w), finally SSIM(A,B)  
%  siga_sq - img A each window gray cov  
%  sigb_sq - img B each window gray cov  
%-------------------------------------------------------------------  
%  Cool_ben  
%========================================================================  
  
w = fspecial('gaussian', 11, 1.5);  %add window
K(1) = 0.01;                      
K(2) = 0.03;                      
L = 2^16-1;       
ima = double(ima);  
imb = double(imb);  
  
C1 = (K(1)*L)^2;  
C2 = (K(2)*L)^2;  
w = w/sum(sum(w));  
  
ua   = filter2(w, ima, 'valid');%Gaussain conv 
ub   = filter2(w, imb, 'valid');  
ua_sq = ua.*ua;  
ub_sq = ub.*ub;  
ua_ub = ua.*ub;  
siga_sq = filter2(w, ima.*ima, 'valid') - ua_sq;  
sigb_sq = filter2(w, imb.*imb, 'valid') - ub_sq;  
sigab = filter2(w, ima.*imb, 'valid') - ua_ub;  
  
ssim_map = ((2*ua_ub + C1).*(2*sigab + C2))./((ua_sq + ub_sq + C1).*(siga_sq + sigb_sq + C2));  
  
  
mssim = mean2(ssim_map);  
  
end  