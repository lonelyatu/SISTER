% This script is the tensor-based dictionary learning using K-CPD
% Key variables:
%           'im'                A multichannel CT image, which can be reconstructed 
%                               using FBP or iterative methods. There is three dimension,
%                              The first and second dimension is the width
%                              and height of the image, and the third
%                              dimension is the spectral dimension.
%
%           'params'        Parameters used in dictionary learning
%
%           'dict'              A four-order tensor. The first three dimension means
%                               the size of an atom. The last dimension means the
%                               number of the atoms.
%
% Yanbo Zhang
% University of Massachusetts Lowell
% yanbozhang007@gmail.com
% 2015-06-26


%% load a multichannel spectral CT image (3D)

clear
clc

% find the file folder of current file
p1 = mfilename('fullpath');
i=strfind (p1,'\');
p1=p1(1:i(end));
cd(p1)

cd('../')                             % one level up from the current folder
codepath = pwd;         
addpath(genpath(codepath))  % add all the sub-folders of this folder to the path

% select image, 1-real 1; 2-real 2; 3-real 3; 4- simulation 1
ind_data = 2; 

% normalized the images so that all channel images have same means
isimaver = 2;   

switch ind_data
    case 1
        % the first real scanned mouse dataset
        load([codepath,'\Results\Mouse_2014-6-1-15-48-5\Alldata'],'ReconrealMC')
        im = ReconrealMC;
    case 2
        % the second real scanned mouse dataset
        load('E:\hudianlin\ÄÜÆ×Îé\Hu_1\FBP\ReconNoise512-180.mat','ReconNoise')
%         ReconNoise = rot90(ReconNoise,3);
%         for i = 1 : size(ReconNoise,3)
%             ReconNoise(:,:,i) = fliplr(ReconNoise(:,:,i));
%         end
        im = ReconNoise;
    case 3
        % the real scanned lamb dataset
        disp('Please select other datasets.')
        return
    case 4
        % the simulated mouse dataset
%         load('D:\Yanbo\code\TDL\DataPrepare\simulation\Alldata_FBP','Recon','ReconNoise')
%         load([codepath,'\DataPrepare\simulation\Alldata_FBP2'],'ReconNoise')
%         load([codepath,'\DataPrepare\simulation\Alldata_FBP_80'],'ReconNoise')
%         load C:\Users\Weiwen_Wu\Desktop\BM4D\4DSR\Data\SART_Simulation_dRcnLamb_0.03_nMulti_2_nPhonExp_4\g_dDcm_Simulation_SART.mat
        load([codepath,'\DataPrepare\simulation\Alldata_FBP_256_256_8_640'],'Recon')
%         im=gather(g_dDcm_SART);
        im = Recon;%Recon;
    otherwise
        error('Please input the right image index !')
end

if isimaver
    % normalization, so that all channels have the same norm
    nChannel = size(im, 3);
    imratio = zeros(1, nChannel);
    immeanall = sum(im(:).^2)/nChannel;
    for i = 1:nChannel
        temp = im(:,:,i);
        immean = sum(temp(:).^2);
        imratio(i) = sqrt(immeanall/immean);      
        im(:,:,i) = im(:,:,i)*imratio(i);
    end
end

%% parameter setting for K-CPD

nChannel = size(im, 3);
params.x = im;
params.imratio = imratio;
params.blocksize = 8;
% number of atoms
params.dictsize = 1024;%1024;512, 1024, 1536 and 2048
% error of sparse coding
params.sigma = 0;%0.20;%0;
% sparsity of dictionary learning, 5-8
params.maxatoms = 5;%5;6
params.maxval = 255;

params.trainnum = 120000;%40000;
params.iternum = 100;
params.memusage = 'high';
params.tensorblocksize = [params.blocksize params.blocksize nChannel];

% iteration number of K-CPD
params.maxiters = 100;
params.tol = 1.0e-06;
params.printitn = 5;
% type of direct currency removement: 'total' removes the total DC of
% a block, 'channel' removes the DCs of each channel in a block
params.dctype = 'channel'; % 'total'£¬'channel'

[params,verbose,msgdelta] = parasetCPD(params);

% select only part of the data that contains strong signal to be the
% training dataset
trainnumSelect = 10000;
params = removenoisedataTensorSelect(params, trainnumSelect);


%% Dictionary initialization

% type of dictionary initialization, 0-DCT; 1-Data; 2-Data clustering
params.inittype = 1;

if params.inittype == 0
    dictdimshow = [10 10 6];% [8 8 4];
    temp1 =odctdict(params.tensorblocksize(1),dictdimshow(1));
    temp2 =odctdict(params.tensorblocksize(2),dictdimshow(2));
    temp3 =odctdict(params.tensorblocksize(3),dictdimshow(3));
    Dinit = zeros([params.tensorblocksize params.dictsize]);
    indxijk = 0;
    
    for dk = 1:dictdimshow(3)
        d3 = temp3(:,dk);
        for dj = 1:dictdimshow(2)
            d2 = temp2(:,dj);
            for di = 1:dictdimshow(1)
                d1 = temp1(:,di);
                indxijk = indxijk+1;
                atomtensor = tensor(ktensor(1,d1, d2,d3));
                atomtemp = atomtensor.data;
                if indxijk <= params.dictsize
                    Dinit(:,:,:,indxijk) = atomtemp/norm(atomtemp(:));
                end
            end
        end
    end
    
elseif params.inittype == 1
    % Initialization method 2: selecting from data directly, although it's not rank-1 at the beginning
    datasize = size(params.data);
    perm = randperm(datasize(end));
    perm = perm(1:params.dictsize);
    Dinit = params.data(:,:,:,perm);    
    
elseif params.inittype == 2
    % Initialization method3: Decomposing the data blocks in rank-n, then clustering them into 512 
    % classes, get the centers of the 512 classes.
    error('This dictionary initialization method is not availible now.')
    
end

% initialized dictionary
params.Dinit = Dinit;   

%% dictionary learning

dict = kcpd20140705(params,verbose,msgdelta);

% figure,displayDictionaryElementsAsImage(dict, 16, 16,8,8,0);
dictemp = squeeze(dict(1:params.tensorblocksize(1)*params.tensorblocksize(2),:));
figure,displayDictionaryElementsAsImage(dictemp,32, params.dictsize/32, params.tensorblocksize(1),params.tensorblocksize(2),0);


%% save the dictionary

% saving path of the learned dictionary
filesavepath0 = [codepath,'\Results\'];
t = clock;
tstr = strcat('_',num2str(t(1)),'-',num2str(t(2)),'-',num2str(t(3)),'-',num2str(t(4)),'-',num2str(t(5)),'-',num2str(floor(t(6))));
% fSaveName = strcat('dict_CPD_ratio_params_case', num2str(ind_data),'_', tstr); 
fSaveName = 'Dict_256_256_8_640';
% save(strcat(filesavepath0,fSaveName),'dict','params'); 
save('E:\hudianlin\ÄÜÆ×Îé\Hu_1\Dict_512_512_8_180_real','dict','params')
