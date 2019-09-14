
reconsize = 256;
params.Dsource2centr = 5;
params.Dsource2detec = 10;
params.NumofView = 720;
params.NumofBin = 512;       % number of detector bins
params.pixelsize = 2/256;
params.binsize = 0.008;
params.binshift = 0;               % detector shift, mm

params.proj = idealSiProj(:,:,1);
iminit = zeros(reconsize);
iternum = 20;
OSnum = 20;

%%

paramsones = params;
paramsones.im = ones(reconsize);
projones = projdd(paramsones);

% backprojection
paramsones.proj = projones;
paramsones.reconsize = reconsize;
bprojones = bprojdd(paramsones);

%%
imupdate = iminit;
figure,
for i = 1:iternum
    tic
    for ios = 1:OSnum
        imupdate = OSSART(params, imupdate, ios, OSnum, bprojones);
        imshow(imupdate,[])
        title(num2str(i))
        drawnow
    end
    toc
end

%% Multichannle OS-SART

nMC = 8;
iminitMC = zeros(reconsize, reconsize, nMC);
imupdateMC = iminitMC;
params.projMC = idealSiProj;
figure,
for i = 1:iternum
    tic
    for ios = 1:OSnum
        imupdateMCcore = OSSARTCoreMCdd(params, imupdateMC, ios, OSnum);
        imupdateMC = imupdateMC +imupdateMCcore./repmat(bprojones,1,1,nMC);
        imshow(imupdateMC(:,:,2),[])
        title(num2str(i))
        drawnow
    end
    toc
end