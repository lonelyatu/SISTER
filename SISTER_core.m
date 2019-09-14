function [imgReconMC] = SISTER_core(ProjMC, projTotalWeightMC, params, Niter_OSnum, CPparams)


reconsize = size(params.imgInitMC,1);
[~, NumofView, nMC] = size(ProjMC);

Niter = Niter_OSnum(1, :);          % vector of main loops numbers
NiterSum = sum(Niter);
OSnum = Niter_OSnum(2, :);      % vectior of subsets numbers
OSnumIter = [];                           % subsets numbers of each iteration
for i = 1:size(Niter_OSnum,2)
    OSnumIter = [OSnumIter, OSnum(i)*ones(1, Niter(i))];
end

lambdaDict = params.lambdaDict;
w=params.w;
k0=params.k0;
cntpartMC = zeros(size(params.imgInitMC));
for iMC = 1:nMC
    cntpartMC(:,:,iMC) = params.cnt;
end

imgkMC = params.imgInitMC;    % initial image
beta=params.beta;
%% compute the denominator of fidelity part

paramsones = params;
paramsones.im = ones(reconsize);
projones = projdd(paramsones);
paramsones.reconsize = reconsize;
imgRaypixsumMC = zeros(reconsize, reconsize, nMC);
for iMC = 1:nMC
    % backprojection
    paramsones.proj = projTotalWeightMC(:, :, iMC).*projones;
    bprojones = bprojdd(paramsones);
    imgRaypixsumMC(:,:,iMC) = bprojones;
end
ratedict = sum(sum(sum(imgRaypixsumMC)))/sum(sum(sum(cntpartMC)));

%% main loop
U = zeros(reconsize, reconsize, nMC);
T=zeros(reconsize, reconsize, nMC);
% normalize each channel
ProjMCratio = imratioMC(ProjMC, params.imratio);
imgkMC = imratioMC(imgkMC, params.imratio);

params.projMC = ProjMCratio;

beta_weights = beta * ones(size(U));
lambda_weights = lambdaDict * ones(size(U));
beta_g = zeros(size(U));
lambdaDict_g = zeros(size(U));


for kIter = 1:NiterSum

    kIter
    OS = OSnumIter(kIter);          % subset number of current iteration
    
    % iteration of each subset
    for iOS = 1:OS
        
        % compute numerator of fidelity part
        indViewOS = iOS:OS:NumofView;
        gradFidMC = - OSSARTCoreMCwdd(params, imgkMC, projTotalWeightMC, indViewOS);
        
        % OS-SART update
        if iOS ~= OS
            imgkMC = imgkMC - NumofView/length(indViewOS)*gradFidMC./imgRaypixsumMC;
        end
        
        %         % nonnegtivity
        imgkMC(imgkMC<0) = 0;
        
    end

    imgkMC0 = imgkMC;
    params.x = imgkMC0;
    imgDict = CP_Cluster(imgkMC0,CPparams);
    weights = computChannelWeight(imgDict);
    for channel = 1 : size(beta_weights,3)
        beta_g(:,:,channel) = beta_weights(:,:,channel) .* weights(channel);
        lambdaDict_g(:,:,channel) = lambda_weights(:,:,channel) .* weights(channel);
    end
    
    imgkMC = imgkMC0 - (OS*gradFidMC + lambdaDict_g.*cntpartMC.*ratedict.*(imgkMC0 - imgDict)+ratedict.*beta_g.*(imgkMC0-U-T))./...
        (ratedict.*beta_g+imgRaypixsumMC + ratedict.*lambdaDict_g.*cntpartMC);

    for channel=1:size(imgkMC,3)          
       U(:,:,channel) =  L0Smoothing(imgkMC(:,:,channel)-T(:,:,channel),w,k0);
       
    end
    T =  T-imgkMC+U;
    loc=find(imgkMC<0);
    imgkMC(loc)=0;

end
imgReconMC = imratioMC(imgkMC, 1./params.imratio);
