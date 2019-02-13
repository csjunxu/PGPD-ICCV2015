%------------------------------------------------------------------------------------------------
% PGPD_Denoising_faster - Denoising by Weighted Sparse Coding 
%                                           with Learned Patch Group Prior.
% CalNonLocal_II - Calculate the non-local similar patches (Noisy Patch Groups)
%                                 by Integral Image (or Summed Area Table, SAT) Techniques
% Author:  Jun Xu, csjunxu@comp.polyu.edu.hk
%              The Hong Kong Polytechnic University
%------------------------------------------------------------------------------------------------
function  [im_out,par] = PGPD_Denoising_faster(par,model)
im_out = par.nim;
par.nSig0 = par.nSig;
par.Win = min(2*par.ps,par.Win); % size of window around the patch
% parameters for noisy image
[h,  w] = size(im_out);
par.maxr = h-par.ps+1;
par.maxc = w-par.ps+1;
par.h = h;
par.w = w;
r = 1:par.step:par.maxr;
par.r = [r r(end)+1:par.maxr];
c = 1:par.step:par.maxc;
par.c = [c c(end)+1:par.maxc];
par.lenr = length(par.r);
par.lenc = length(par.c);
par.ps2 = par.ps^2;
par.maxrc = par.maxr*par.maxc;
par.lenrc = par.lenr*par.lenc;
% parameters for Pad noisy image
hp = h + 2*par.Win;
wp = w + 2*par.Win;
par.maxrp = hp-par.ps+1;
par.maxcp = wp-par.ps+1;
par.maxrcp = par.maxrp*par.maxcp;
% positions for integral image
par.x = par.r+par.ps+par.Win-1;
par.y = par.c+par.ps+par.Win-1;
par.x0 = par.r+par.Win-1;
par.y0 = par.c+par.Win-1;
for ite = 1 : par.IteNum
    % iterative regularization
    im_out = im_out+par.delta*(par.nim-im_out);
    % estimation of noise variance
    if ite == 1
        par.nSig = par.nSig0;
    else
        dif = mean( mean( (par.nim-im_out).^2 ) ) ;
        par.nSig = sqrt( abs( par.nSig0^2 - dif ) )*par.eta;
    end
    % search non-local patch groups
    [nDCnlX,blk_arr,DC,par] = Image2PG_II( im_out, par);
    % Gaussian dictionary selection by MAP
    if mod(ite-1,2) == 0
        PYZ = zeros(model.nmodels,size(DC,2));
        sigma2I = par.nSig^2*eye(par.ps2);
        for i = 1:model.nmodels
            sigma = model.covs(:,:,i) + sigma2I;
            [R,~] = chol(sigma);
            Q = R'\nDCnlX;
            TempPYZ = - sum(log(diag(R))) - dot(Q,Q,1)/2;
            TempPYZ = reshape(TempPYZ,[par.nlsp size(DC,2)]);
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
    X_hat = zeros(par.ps2,par.maxrc,'single');
    W = zeros(par.ps2,par.maxrc,'single');
    for   j = 1:length(seg)-1
        idx = s_idx(seg(j)+1:seg(j+1));
        cls = dicidx(idx(1));
        D   = par.D(:,:, cls);
        S    = par.S(:,cls);
        lambdaM = repmat(par.c1*par.nSig^2./ (sqrt(S)+eps),[1 par.nlsp]);
        for i = 1:size(idx,1)
            Y = nDCnlX(:,(idx(i)-1)*par.nlsp + 1:idx(i)*par.nlsp);
            b = D'*Y;
            % soft threshold
            alpha = sign(b).*max(abs(b) - lambdaM/2,0);
            % add DC components and aggregation
            X_hat(:,blk_arr(:,idx(i))) = X_hat(:,blk_arr(:,idx(i))) + bsxfun(@plus,D*alpha, DC(:,idx(i)));
            W(:,blk_arr(:,idx(i))) = W(:,blk_arr(:,idx(i))) + ones(par.ps2,par.nlsp);
        end
    end
    % Reconstruction
    im_out = zeros(h,w,'single');
    im_wei = zeros(h,w,'single');
    r = 1:par.maxr;
    c = 1:par.maxc;
    k = 0;
    for i = 1:par.ps
        for j = 1:par.ps
            k = k + 1;
            im_out(r-1+i,c-1+j) = im_out(r-1+i,c-1+j) + reshape( X_hat(k,:)', [par.maxr par.maxc]);
            im_wei(r-1+i,c-1+j) = im_wei(r-1+i,c-1+j) + reshape( W(k,:)', [par.maxr par.maxc]);
        end
    end
    im_out = im_out./im_wei;
    % calculate the PSNR and SSIM
    PSNR = csnr( im_out*255, par.I*255, 0, 0 );
    SSIM = cal_ssim( im_out*255, par.I*255, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n',ite, PSNR,SSIM);
end
im_out(im_out > 1) = 1;
im_out(im_out < 0) = 0;
return; 
