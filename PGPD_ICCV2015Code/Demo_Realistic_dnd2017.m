%--------------------------------------------------------------------------
clear;
addpath('model');
addpath('NoiseEstimation');
Original_image_dir = 'dnd_2017/images_srgb/';
fpath = fullfile(Original_image_dir, '*.mat');
im_dir  = dir(fpath);
im_num = length(im_dir);
load 'dnd_2017/info.mat';

method = 'PGPD';
dataset = 'dnd_2017';
write_MAT_dir = [dataset '_Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

RunTime = [];
for i = 1 :im_num
    load(fullfile(Original_image_dir, im_dir(i).name));
    S = regexp(im_dir(i).name, '\.', 'split');
    [h,w,ch] = size(InoisySRGB);
    for j = 1:size(info(1).boundingboxes,1)
        IMinname = [S{1} '_' num2str(j)];
        bb = info(i).boundingboxes(j,:);
        nim = InoisySRGB(bb(1):bb(3), bb(2):bb(4),:);
        fprintf('%s: \n', IMinname);
        % noise estimation
        t1=clock;
        IMout = zeros(size(nim));
        for c = 1:ch
            nSig = NoiseEstimation(nim(:, :, c)*255, 8);
            % denoising
            [par, model]  =  Parameters_Setting( nSig );
            par.I = nim(:,:,c);
            par.nim = nim(:,:,c);
            [IMoutc,par]  =  PGPD_Denoising_Real(par,model);
            IMout(:,:,c) = IMoutc;
        end
        t2=clock;
        etime(t2,t1)
        alltime(i)  = etime(t2, t1);
        %% output
        IMoutname = sprintf([write_sRGB_dir '/' method '_DND_' IMinname '.png']);
        imwrite(IMout, IMoutname);
    end
end
