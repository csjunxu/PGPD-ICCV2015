%--------------------------------------------------------------------------
clear;
addpath('model');
addpath('NoiseEstimation');
GT_Original_image_dir = 'Real_ccnoise_denoised_part/';
GT_fpath = fullfile(GT_Original_image_dir, '*mean.png');
TT_Original_image_dir = 'Real_ccnoise_denoised_part/';
TT_fpath = fullfile(TT_Original_image_dir, '*real.png');


GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

method = 'PGPD';
dataset = 'dnd2017';
write_MAT_dir = [dataset '_Results/'];
write_sRGB_dir = [write_MAT_dir method];
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
    fprintf('%s \nThe initial PSNR = %2.4f, SSIM = %2.4f. \n', IMname, csnr(IM_GT*255, IM*255, 0, 0 ), cal_ssim(IM_GT*255, IM*255, 0, 0 ));
    for cc = 1:ch
        %% noise estimation
        nSig = NoiseEstimation(IM(:, :, cc)*255, 8);
        %% set parameters
        [par, model]  =  Parameters_Setting( nSig );
        par.I = IM_GT(:,:,cc);
        par.nim = IM(:,:,cc);
        %% denoising
        [IMoutcc,par]  =  PGPD_Denoising_Real(par,model);
        IMout(:,:,cc) = IMoutcc;
    end
    RunTime = [RunTime etime(clock,time0)];
    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
    PSNR = [PSNR csnr( IMout*255, IM_GT*255, 0, 0 )];
    SSIM = [SSIM cal_ssim( IMout*255, IM_GT*255, 0, 0 )];
    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', PSNR(end), SSIM(end));
    imwrite(IMout, [write_sRGB_dir '/' method '_' dataset '_' IMname '.png']);
end
mPSNR = mean(PSNR);
mSSIM = mean(SSIM);
mRunTime = mean(RunTime);
matname = sprintf([write_MAT_dir method '_' dataset '.mat']);
save(matname,'PSNR','mPSNR','SSIM','mSSIM','RunTime','mRunTime');
