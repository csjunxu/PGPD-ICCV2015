%------------------------------------------------------------------------------------------------
% PGPD_Denoising - Denoising by Weighted Sparse Coding
%                                with Learned Patch Group Prior.
% CalNonLocal - Calculate the non-local similar patches (Noisy Patch Groups)
% Author:  Jun Xu, csjunxu@comp.polyu.edu.hk
%              The Hong Kong Polytechnic University
%------------------------------------------------------------------------------------------------
function  [im_out,par] = PGPD_Denoising(par,model)
im_out = par.nim;
par.nSig0 = par.nSig;
[h, w, ch] = size(im_out);
par.maxr = h-par.ps+1;
par.maxc = w-par.ps+1;
par.maxrc = par.maxr * par.maxc;
par.h = h;
par.w = w;
par.ch = ch;
r = 1:par.step:par.maxr;
par.r = [r r(end)+1:par.maxr];
c = 1:par.step:par.maxc;
par.c = [c c(end)+1:par.maxc];
par.lenr = length(par.r);
par.lenc = length(par.c);
par.lenrc = par.lenr*par.lenc;
par.ps2 = par.ps^2;
par.ps2ch = par.ps2*par.ch;
for ite = 1 : par.IteNum
    % iterative regularization
    im_out = im_out + par.delta*(par.nim - im_out);
    % estimate noise variance of each channel
    sigma2I = zeros(par.ps2ch,1);
    for c = 1:ch
        if ite == 1
            par.nSig(c) = par.nSig0(c);
        else
            dif = mean( mean( (par.nim(:,:,c) - im_out(:,:,c)).^2 ) ) ;
            par.nSig(c) = sqrt( abs( par.nSig0(c)^2 - dif ) )*par.eta(c);
        end
        sigma2I((c-1)*par.ps2+1:c*par.ps2,1) = par.nSig(c)^2*ones(par.ps2,1);
    end
    % search non-local patch groups
    [nDCnlX,blk_arr,DC,par] = CalNonLocal( im_out, par);
    % Gaussian dictionary selection by MAP
    if mod(ite-1,2) == 0
        PYZ = zeros(model.nmodels,size(DC,2)/par.nlsp);
        for i = 1:model.nmodels
            sigma = model.covs(:,:,i) + diag(sigma2I);
            [R,~] = chol(sigma);
            Q = R'\nDCnlX;
            TempPYZ = - sum(log(diag(R))) - dot(Q,Q,1)/2;
            TempPYZ = reshape(TempPYZ,[par.nlsp size(DC,2)/par.nlsp]);
            PYZ(i,:) = sum(TempPYZ);
        end
        % find the most likely component for each patch group
        [~,dicidx] = max(PYZ);
        dicidx = dicidx';
        [idx,  s_idx] = sort(dicidx);
        idx2 = idx(1:end-1) - idx(2:end);
        seq = find(idx2);
        seg = [0; seq; length(dicidx)];
    end
    % Weighted Sparse Coding
    X_hat = zeros(par.ps2ch,par.maxrc,'single');
    W = zeros(par.ps2ch,par.maxrc,'single');
    for j = 1:length(seg)-1
        idx =   s_idx(seg(j)+1:seg(j+1));
        idxs = [];
        for i = 1:length(idx)
            idxs = [idxs (idx(i)-1)*par.nlsp+1:idx(i)*par.nlsp];
        end
        cls =   dicidx(idx(1));
        D   =   par.D(:,:, cls);
        S    = par.S(:,cls);
        lambdaM = repmat(par.c1*sigma2I./ (sqrt(S)+eps ),[1 length(idxs)]);
        Y = nDCnlX(:,idxs);
        b = D'*Y;
        % soft threshold
        alpha = sign(b).*max(abs(b)-lambdaM/2,0);
        % add DC components and aggregation
        X_hat(:,blk_arr(:,idxs)) = X_hat(:,blk_arr(:,idxs))+bsxfun(@plus,D*alpha, DC(:,idxs));
        W(:,blk_arr(:,idxs)) = W(:,blk_arr(:,idxs))+ones(par.ps2ch, length(idxs));
    end
    % Reconstruction
    im_out = zeros(h,w,ch,'single');
    im_wei = zeros(h,w,ch,'single');
    r = 1:par.maxr;
    c = 1:par.maxc;
    k = 0;
    for l = 1:1:par.ch
        for i = 1:par.ps
            for j = 1:par.ps
                k = k+1;
                im_out(r-1+i, c-1+j, l)  =  im_out(r-1+i, c-1+j, l) + reshape( X_hat(k,:)', [par.maxr par.maxc] );
                im_wei(r-1+i, c-1+j, l)  =  im_wei(r-1+i, c-1+j, l) + reshape( W(k,:)', [par.maxr par.maxc] );
            end
        end
    end
    im_out  =  im_out./im_wei;
    % calculate the PSNR and SSIM
    PSNR = csnr( im_out*255, par.I*255, 0, 0 );
    SSIM = cal_ssim( im_out*255, par.I*255, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n',ite, PSNR,SSIM);
end
im_out(im_out > 1) = 1;
im_out(im_out < 0) = 0;
return;