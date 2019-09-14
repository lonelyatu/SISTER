% projection and backprojection based on distance-driven method
%
% 2015-10-09


% clear all
% close all
% clc

im = phantom(256);%phantom(256);%ones(256);
reconsize = size(im,1);
NumofView = 720;
startAngle = 0;
endAngle = 360-360/NumofView;
IsEquiDis = 1;          % 1: equidistance detector; 0: equi-angel detector
Dsource2centr = 5;
Dsource2detec = 10;
binsize = 0.008;
pixelsize = 2/256;%1.5;

binshift = 0;               % detector shift, mm
NumofBin = 512;       % number of detector bins

%% generate projection

params.im = im;
params.Dsource2centr = Dsource2centr;
params.Dsource2detec = Dsource2detec;
params.NumofView = NumofView;
params.NumofBin = NumofBin;       % number of detector bins
params.pixelsize = pixelsize;
params.binsize = binsize;
params.binshift = binshift;               % detector shift, mm

tic
proj = projdd(params);
toc
figure,imshow(proj,[])

%% backprojection

clear params
params.proj = proj;
params.reconsize = reconsize;
params.Dsource2centr = Dsource2centr;
params.Dsource2detec = Dsource2detec;
[NumofBin, NumofView] = size(proj);
params.pixelsize = pixelsize;
params.binsize = binsize;
params.binshift = binshift; 

tic
bim = bprojdd(params);
toc
figure,imshow(bim,[])

%%  generate ones projection and backprojection

% reconsize = 256;
params.im = ones(reconsize);
params.Dsource2centr = Dsource2centr;
params.Dsource2detec = Dsource2detec;
params.NumofView = NumofView;
params.NumofBin = NumofBin;       % number of detector bins
params.pixelsize = pixelsize;
params.binsize = binsize;
params.binshift = binshift;               % detector shift, mm

tic
% projection
projones = projdd(params);

% backprojection
params.proj = projones;
params.reconsize = reconsize;
bprojones = bprojdd(params);
toc

figure,imshow(bprojones,[])

%%

bimnorm = bim./bprojones;
params.im = bimnorm;
proj2 = projdd(params);
figure,imshow(proj2,[])

%% SART

iternum = 20;
imupdate = zeros(reconsize);
figure,
for i = 1:iternum
    
    params.im = imupdate;
    proj2 = projdd(params);
    
    params.proj = proj - proj2;
    imdiff = bprojdd(params);
    imupdate = imupdate + imdiff./bprojones;
    imshow(imupdate,[])
    title(num2str(i))
    drawnow 
end


%% OS-SART

iternum = 20;
OSnum = 20;
imupdate = zeros(reconsize);
params.iViews = 0 : 360/NumofView : 360 - 360/NumofView;
figure,
for i = 1:iternum 
    tic
    for ios = 1:OSnum
        indViewOS = ios:OSnum:NumofView;
        params.iViews = params.iViews(indViewOS);
%         params.iViews = indViewOS;
        params.im = imupdate;
        proj2OS = projdd(params);
        
        params.proj = proj(:,indViewOS) - proj2OS;
        imdiff = bprojdd(params);
        imupdate = imupdate + NumofView/length(indViewOS)*imdiff./bprojones;
    end
    toc
    imshow(imupdate,[])
    title(num2str(i))
    drawnow
end
%% FBP 


addpath('D:\Yanbo\code\Spectral CT real data processing\')
NumofRou = NumofBin;%1024;
dectsize = binsize;%  计算得到的一个detector bin的尺寸为0.0553mm，理论值为0.0550mm
Dsourc2obj = Dsource2centr;
Dsourc2detec = Dsource2detec; 
Dobj2detec =  Dsourc2detec - Dsourc2obj;

Mirror = zeros(reconsize,reconsize);
NumofTheta = NumofView;%660;
NumofRoudown = round(NumofRou);  
dudown = dectsize;
DofSourceCentr = Dsourc2obj;
DofSourceDetector = Dsourc2detec;

dectectorshift = 1;                     % 该值为1表示检测器单元左右顺序颠倒，否则为0

YMin = -1;                              % 此4值为重建图像范围
YMax = 1;
XMin =  -1;
XMax = 1;

MagFactor = 1;

tic
Proj1 = reshape(flipud(proj),NumofRoudown,1,NumofView);
Recon = Recon_offcenter_fullscan_pixsize_MTF_rot(Proj1,YMin,YMax,XMin,XMax,DofSourceCentr,Dobj2detec,dudown,dectsize,pixelsize,reconsize,0,360,dectectorshift);
toc
figure,imshow(Recon,[])

