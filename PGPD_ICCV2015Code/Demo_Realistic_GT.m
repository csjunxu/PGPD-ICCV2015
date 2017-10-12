%--------------------------------------------------------------------------
clear;
addpath('model');
addpath('NoiseEstimation');
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\DJI_Results\Real_MeanImage\';
% GT_fpath = fullfile(GT_Original_image_dir, '*.JPG');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\DJI_Results\Real_NoisyImage\';
% TT_fpath = fullfile(TT_Original_image_dir, '*.JPG');
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_MeanImage\';
% GT_fpath = fullfile(GT_Original_image_dir, '*.png');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_NoisyImage\';
% TT_fpath = fullfile(TT_Original_image_dir, '*.png');
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_ccnoise_denoised_part\';
% GT_fpath = fullfile(GT_Original_image_dir, '*mean.png');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_ccnoise_denoised_part\';
% TT_fpath = fullfile(TT_Original_image_dir, '*real.png');
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\our_Results\Real_MeanImage\';
% GT_fpath = fullfile(GT_Original_image_dir, '*.JPG');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\our_Results\Real_NoisyImage\';
% TT_fpath = fullfile(TT_Original_image_dir, '*.JPG');
GT_Original_image_dir = 'C:/Users/csjunxu/Desktop/RID_Dataset/RealisticImage/';
GT_fpath = fullfile(GT_Original_image_dir, '*mean.JPG');
TT_Original_image_dir = 'C:/Users/csjunxu/Desktop/RID_Dataset/RealisticImage/';
TT_fpath = fullfile(TT_Original_image_dir, '*real.JPG');

GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

method = 'PGPD';
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/PolyU_Results/'];
write_sRGB_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/PolyU_Results/' method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end
PSNR = [];
SSIM = [];
RunTime = [];
for i = 1 : im_num
    IM =   im2double(imread( fullfile(TT_Original_image_dir,TT_im_dir(i).name) ));
    IM_GT = im2double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
    % S = regexp(TT_im_dir(i).name, '\.', 'split');
    IMname = TT_im_dir(i).name(1:end-9);
    [h,w,ch] = size(IM);
    time0 = clock;
    IMout = zeros(size(IM));
    for cc = 1:ch
        %% noise estimation
        nSig = NoiseEstimation(IM(:, :, cc)*255, 8);
        %% set parameters
        [par, model]  =  Parameters_Setting( nSig );
        par.I = IM_GT(:,:,cc);
        par.nim = IM(:,:,cc);
        fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', csnr(IM_GT(:,:,cc)*255, IM(:,:,cc)*255, 0, 0 ), cal_ssim(IM_GT(:,:,cc)*255, IM(:,:,cc)*255, 0, 0 ));
        %% denoising
        [IMoutcc,par]  =  PGPD_Denoising_Real(par,model);
        IMout(:,:,cc) = IMoutcc;
    end
    RunTime = [RunTime etime(clock,time0)];
    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
    PSNR = [PSNR csnr( IMout*255, IM_GT*255, 0, 0 )];
    SSIM = [SSIM cal_ssim( IMout*255, IM_GT*255, 0, 0 )];
    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', PSNR(end), SSIM(end));
    imwrite(IMout, [write_sRGB_dir '/' method '_our_' IMname '.png']);
end
mPSNR = mean(PSNR);
mSSIM = mean(SSIM);
mRunTime = mean(RunTime);
matname = sprintf([write_MAT_dir method '_our.mat']);
save(matname,'PSNR','mPSNR','SSIM','mSSIM','RunTime','mRunTime');