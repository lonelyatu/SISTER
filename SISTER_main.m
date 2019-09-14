
%% addpath

clear
clc

p1 = mfilename('fullpath');
i=strfind (p1,'\');
p1=p1(1:i(end));
cd(p1)

% cd('../')                             % one level up from the current folder
codepath = pwd;         
addpath(genpath(codepath))  % add all the sub-folders of this folder to the path


% the simulated mouse dataset
load('./demo_data','ProjEnergyNoise','ReconNoise', 'pixelsize', 'cnt')
im = ReconNoise;
ProjMC = zeros(size(ProjEnergyNoise));
nChannel = size(ProjEnergyNoise, 3);

for i = 1:nChannel
    ProjMC(:,1,i) = ProjEnergyNoise(:,1,i);
    ProjMC(:,2:end,i) = fliplr(ProjEnergyNoise(:,2:end,i));
end

imratio = zeros(1, nChannel);
immeanall = sum(im(:).^2)/nChannel;
for i = 1:nChannel
    temp = im(:,:,i);
    immean = sum(temp(:).^2);
    imratio(i) = sqrt(immeanall/immean);
end

params.imratio = imratio;
params.projMC = ProjMC;
params.Dsource2centr = 132;
params.Dsource2detec = 180;
params.NumofView = 160;
params.NumofBin = 512;       % number of detector bins
params.pixelsize = 0.075*2;
params.binsize = 0.1;
params.binshift = 0;               % detector shift, mm


projBadDetecWeightMC = ones(size(ProjMC)); 
projTotalWeightMC = ones(size(ProjMC)); 

params.reconsize = 256;
imgInitMC = zeros(params.reconsize, params.reconsize, nChannel);   % initalized image

% Niter_OSnum is a two rows matrix, which specifies the numbers of main loops 
% and corresponding subsets respectively of OS
Niter_OSnum = [100; 10];

params.imgInitMC = imgInitMC;  
params.memusage = 'high'; 

params.cnt = cnt;
params.k0=1.1;
params.lambdaDict = 1.1;
params.w = 1.8*1e-4;
params.beta = 5.7;

SISTERparams.atoms = 32;
SISTERparams.error = 1e-4;
SISTERparams.k = 50;
SISTERparams.class = 100;

Recon_SISTER = SISTER_core(ProjMC, projTotalWeightMC, params, Niter_OSnum, SISTERparams);


