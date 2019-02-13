function  [im_out,par] = PGPD_Denoising_Real(par,model)
im_out = par.nim;
par.nSig0 = par.nSig;
[h,  w, ch] = size(im_out);
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
fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', csnr(par.nim*255, par.I*255, 0, 0 ), cal_ssim(par.nim*255, par.I*255, 0, 0 ));
for ite = 1 : par.IteNum
    % iterative regularization
    im_out = im_out + par.delta*(par.nim - im_out);
    % estimation of noise variance
    if ite == 1
        par.nSig = par.nSig0;
    else
        dif = mean( mean( (par.nim - im_out).^2 ) ) ;
        par.nSig = sqrt( abs( par.nSig0^2 - dif ) )*par.eta;
    end
    % search non-local patch groups
    [nDCnlX,blk_arr,DC,par] = Image2PG( im_out, par);
    % Gaussian dictionary selection by MAP
        if mod(ite-1,2) == 0
        %% GMM: full posterior calculation
        nPG = size(nDCnlX,2)/par.nlsp; % number of PGs
        PYZ = zeros(model.nmodels,nPG);
        for i = 1:model.nmodels
            sigma = model.covs(:,:,i);
            [R,~] = chol(sigma);
            Q = R'\nDCnlX;
            TempPYZ = - sum(log(diag(R))) - dot(Q,Q,1)/2;
            TempPYZ = reshape(TempPYZ,[par.nlsp nPG]);
            PYZ(i,:) = sum(TempPYZ);
        end
        %% find the most likely component for each patch group
        [~,dicidx] = max(PYZ);
        dicidx=repmat(dicidx, [par.nlsp 1]);
        dicidx = dicidx(:);
        [idx,  s_idx] = sort(dicidx);
        idx2 = idx(1:end-1) - idx(2:end);
        seq = find(idx2);
        seg = [0; seq; length(dicidx)];
        end
    
    % Weighted Sparse Coding
    X_hat = zeros(par.ps2ch,par.maxrc,'double');
    W_hat = zeros(par.ps2ch,par.maxrc,'double');
    for   j = 1:length(seg)-1
        idx =   s_idx(seg(j)+1:seg(j+1));
        cls =   dicidx(idx(1));
        D   =   par.D(:,:, cls);
        S    = par.S(:,cls);
        Y = nDCnlX(:,idx);
        lambdaM = repmat(par.c1*par.nSig^2./ (sqrt(S)+eps ),[1 length(idx)]);
        b = D'*Y;
        % soft threshold
        alpha = sign(b).*max(abs(b)-lambdaM/2,0);
        % add DC components and aggregation
        X_hat(:,blk_arr(:,idx)) = X_hat(:,blk_arr(:,idx))+bsxfun(@plus,D*alpha, DC(:,idx));
        W_hat(:,blk_arr(:,idx)) = W_hat(:,blk_arr(:,idx))+ones(par.ps2ch, length(idx));
    end
    % Reconstruction
    im_out = PGs2Image(X_hat,W_hat,par);
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out*255, par.I*255, 0, 0 );
    SSIM      =  cal_ssim( im_out*255, par.I*255, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n',ite, PSNR,SSIM);
end
im_out(im_out > 1) = 1;
im_out(im_out < 0) = 0;
return;
