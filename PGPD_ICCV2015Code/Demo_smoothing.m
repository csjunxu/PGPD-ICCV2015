clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\L0smoothing\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, 'pflower.jpg');
im_dir  = dir(fpath);
im_num = length(im_dir);

method = 'PGPD';
load './model/PG_GMM_8x8_win15_nlsp10_delta0.002_cls33.mat';
nSig = 10;
par.step = 3;       % the step of two neighbor patches
par.IteNum = 4;  % the iteration number
par.nSig      =   nSig/255;
par.ps = ps;        % patch size
par.nlsp = nlsp;  % number of non-local patches
par.Win = win;   % size of window around the patch
% dictionary and regularization parameter
for i = 1:size(GMM_D,2)
    par.D(:,:,i) = reshape(single(GMM_D(:, i)), size(GMM_S,1), size(GMM_S,1));
end
par.S = single(GMM_S);
writefilepath = 'C:/Users/csjunxu/Desktop/L0smoothing/';

delta = 0;
par.delta = delta;
for c1  = .1:.1:1
    par.c1 = c1*2*sqrt(2);
    for eta = 1
        par.eta = eta;
        for i = 1:im_num
            % set parameters
            [par, model]  =  Parameters_Setting( nSig );
            % read image
            nim = im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            nim_ycbcr = rgb2ycbcr(nim);
            nim_y = nim_ycbcr(:,:,1);
            nim_cb = nim_ycbcr(:,:,2);
            nim_cr = nim_ycbcr(:,:,3);
            par.nim = nim_y;
            par.I = par.nim;
            % PGPD denoising
            [rim_y,par]  =  PGPD_Denoising_faster(par,model); % faster speed
            nim_ycbcr(:,:,1) = rim_y;
            rim = ycbcr2rgb(nim_ycbcr);
            imname = sprintf([writefilepath 'pflower_' method '_nSig' num2str(nSig)  '_delta' num2str(delta) '_c1_' num2str(c1) '_eta' num2str(eta) '.png']);
            imwrite(rim,imname);
        end
    end
end
