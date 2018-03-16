clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\TWSCGIN\cleanimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

method = 'PGPD';
writematpath = 'C:/Users/csjunxu/Desktop/ThesisComments/Results_AWGN/';
writefilepath  = [writematpath method '/'];
if ~isdir(writefilepath)
    mkdir(writefilepath);
end
load './model/PG_GMM_9x9_win15_nlsp10_delta0.002_cls33.mat';
for i = 1:size(GMM_D,2)
    par.D(:,:,i) = reshape(single(GMM_D(:, i)), size(GMM_S,1), size(GMM_S,1));
end
par.S = single(GMM_S);
for nSig = [150 200]
    for delta = [0.03]
        par.delta = delta;
        for c1 = .05
            par.c1 = c1*2*sqrt(2);
            for eta = .8:.1:1.5
                par.eta=eta;
                PSNR = [];
                SSIM = [];
                for i = 1:im_num
                    % set parameters
                    par.step = 3;       % the step of two neighbor patches
                    par.IteNum = 4;  % the iteration number
                    par.nSig      =   nSig/255;
                    par.ps = ps;        % patch size
                    par.nlsp = nlsp;  % number of non-local patches
                    par.Win = win;   % size of window around the patch
                    % read clean image
                    par.I = im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                    % generate noisy image
                    randn('seed',0);
                    par.nim =   par.I + par.nSig*randn(size(par.I));
                    fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', csnr( par.nim*255, par.I*255, 0, 0 ),cal_ssim( par.nim*255, par.I*255, 0, 0 ));
                    % PGPD denoising
                    [im_out,par]  =  PGPD_Denoising(par,model);
                    % [im_out,par]  =  PGPD_Denoising_faster(par,model); % faster speed
                    %                     imname = sprintf([writefilepath method '_nSig' num2str(nSig)  '_delta' num2str(delta) '_c1_' num2str(c1) '_eta' num2str(eta) '_delta' num2str(delta) '_' im_dir(i).name]);
                    %                     imwrite(im_out,imname);
                    % calculate the PSNR and SSIM
                    PSNR = [PSNR csnr( im_out*255, par.I*255, 0, 0 )];
                    SSIM =  [SSIM cal_ssim( im_out*255, par.I*255, 0, 0 )];
                    fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n', im_dir(i).name, PSNR(end), SSIM(end));
                end
                mPSNR=mean(PSNR);
                mSSIM=mean(SSIM);
                name = sprintf([writematpath method '_nSig' num2str(nSig) '_delta' num2str(delta) '_c1_' num2str(c1) '_eta' num2str(eta) '.mat']);
                save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM');
            end
        end
    end
end


method = 'PPD';
load './model/PG_GMM_9x9_win15_nlsp10_delta0.002_cls33.mat';
for i = 1:size(GMM_D,2)
    par.D(:,:,i) = reshape(single(GMM_D(:, i)), size(GMM_S,1), size(GMM_S,1));
end
par.S = single(GMM_S);
for nSig = [150 200]
    for delta = [0.03]
        par.delta = delta;
        for c1 = .05
            par.c1 = c1*2*sqrt(2);
            for eta = .8:.1:1.5
                par.eta=eta;
                PSNR = [];
                SSIM = [];
                for i = 1:im_num
                    % set parameters
                    par.step = 3;       % the step of two neighbor patches
                    par.IteNum = 4;  % the iteration number
                    par.nSig      =   nSig/255;
                    par.ps = ps;        % patch size
                    par.nlsp = nlsp;  % number of non-local patches
                    par.Win = win;   % size of window around the patch
                    par.nlsp = 1;
                    % read clean image
                    par.I = im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                    % generate noisy image
                    randn('seed',0);
                    par.nim =   par.I + par.nSig*randn(size(par.I));
                    fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', csnr( par.nim*255, par.I*255, 0, 0 ),cal_ssim( par.nim*255, par.I*255, 0, 0 ));
                    % PGPD denoising
                    [im_out,par]  =  PGPD_Denoising(par,model);
                    % [im_out,par]  =  PGPD_Denoising_faster(par,model); % faster speed
                    %                     imname = sprintf([writefilepath method '_nSig' num2str(nSig)  '_delta' num2str(delta) '_c1_' num2str(c1) '_eta' num2str(eta) '_delta' num2str(delta) '_' im_dir(i).name]);
                    %                     imwrite(im_out,imname);
                    % calculate the PSNR and SSIM
                    PSNR = [PSNR csnr( im_out*255, par.I*255, 0, 0 )];
                    SSIM =  [SSIM cal_ssim( im_out*255, par.I*255, 0, 0 )];
                    fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n', im_dir(i).name, PSNR(end), SSIM(end));
                end
                mPSNR=mean(PSNR);
                mSSIM=mean(SSIM);
                name = sprintf([writematpath method '_nSig' num2str(nSig) '_delta' num2str(delta) '_c1_' num2str(c1) '_eta' num2str(eta) '.mat']);
                save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM');
            end
        end
    end
end