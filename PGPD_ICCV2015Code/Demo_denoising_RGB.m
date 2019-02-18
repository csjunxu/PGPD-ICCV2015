%--------------------------------------------------------------------------------------------------
% This is an implementation of the PGPD algorithm for image denoising.
% Author:  Jun Xu, csjunxu@comp.polyu.edu.hk
%              The Hong Kong Polytechnic University
% Please refer to the following paper if you use this code:
% Jun Xu, Lei Zhang, Wangmeng Zuo, David Zhang, and Xiangchu Feng,
% Patch Group Based Nonlocal Self-Similarity Prior Learning for Image Denoising.
% IEEE Int. Conf. Computer Vision (ICCV), Santiago, Chile, December 2015.
% Please see the file License.txt for the license governing this code.
%--------------------------------------------------------------------------------------------------
clear;
Original_image_dir  =    '/home/csjunxu/Github/PGPD_Results/kodak24/';
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

method = 'PGPD';
dataset = 'kodak24';
write_MAT_dir = ['/home/csjunxu/Github/PGPD_Results/' dataset '_Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

nSig = [40 20 30];

% fixed parameters
par.step   = 3;    % the step of two neighbor patches
par.IteNum = 4;  % the iteration number
par.nSig   = nSig/255;
load '/home/csjunxu/Github/PGPD_Results/PGGMM Models/RGB/PG_GMM_RGB_8x8_win30_nlsp10_delta0.002_cls513.mat';
par.ps = ps;        % patch size
par.nlsp = nlsp;  % number of non-local patches
par.Win = win;   % size of window around the patch
% dictionary and regularization parameter
for i = 1:size(GMM_D,2)
    par.D(:,:,i) = reshape(single(GMM_D(:, i)), size(GMM_S,1), size(GMM_S,1));
end
par.S = single(GMM_S);

% tunable parameters
for delta = 0
    par.delta = delta;
    for c1 = .8:.2:2
        par.c1 = c1;
        for etaR = [1 .9 1.1]
            for etaG = [.7 .8 .6]
                for etaB = [.9 1 1.1]
                    par.eta=[etaR etaG etaB];
                    % record all the results in each iteration
                    PSNR = zeros(1, im_num, 'single');
                    SSIM = zeros(1, im_num, 'single');
                    for i = 1:im_num
                        par.nSig0 = nSig/255;
                        % read clean image
                        par.I =  im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                        S = regexp(im_dir(i).name, '\.', 'split');
                        [h, w, ch] = size(par.I);
                        % generate noisy image
                        par.nim = zeros(size(par.I));
                        for c = 1:ch
                            randn('seed',0);
                            par.nim(:, :, c) = par.I(:, :, c) + par.nSig0(c) * randn(size(par.I(:, :, c)));
                        end
                        fprintf('Initial PSNR = %2.4f, SSIM = %2.4f \n', csnr( par.nim*255, par.I*255, 0, 0 ),cal_ssim( par.nim*255, par.I*255, 0, 0 ));
                        % PGPD denoising
                        [im_out,par]  =  PGPD_Denoising(par,model);
                        % calculate the PSNR and SSIM
                        PSNR(i) = csnr( im_out*255, par.I*255, 0, 0 );
                        SSIM(i) = cal_ssim( im_out*255, par.I*255, 0, 0 );
                        fprintf('Denoised PSNR = %2.4f, SSIM = %2.4f\n',PSNR(i),SSIM(i));
                    end
                    mPSNR=mean(PSNR);
                    mSSIM=mean(SSIM);
                    matname = sprintf([write_MAT_dir method '_' num2str(cls_num) '_' dataset '_' ...
                        num2str(c1) '_' num2str(etaR) '_' num2str(etaG) '_' num2str(etaB) '.mat']);
                    save(matname,'PSNR','SSIM','mPSNR','mSSIM');
                end
            end
        end
    end
end