function imupdateMC = OSSARTCoreMCwdd(params, iminitMC, projTotalWeightMC, indViewOS)
%
% The core of weighted OSSART for multichannel
%
% Input variables:
%           'params'                            structure dataset
%           'projTotalWeightMC'       projection weighting, 3D
%           'iminitMC'                         input image, multichannel, 3D
%           'indViewOS'                      current used index of projection view
%
% Output variables:
%           'imupdateMC'                   updated result
%
% Yanbo Zhang
% University of Massachusetts Lowell
% yanbozhang007@gmail.com
% 2015-10-14

imupdateMC = iminitMC; %  现在的迭代结果图像
projMC = params.projMC; % 投影数据
nMC = size(projMC,3);   % 计算现在的重建图像投影数据和标准投影数据之间的差异，然后然反投影回来 A_T(Ax_s - y_s)

iViews = 0 : 360/params.NumofView : 360 - 360/params.NumofView;
for iMC = 1:nMC
    
    params.iViews = iViews(indViewOS);
    params.im = imupdateMC(:,:,iMC);
    proj2OS = projdd(params);
    
    params.proj = projTotalWeightMC(:,indViewOS,iMC).*(projMC(:,indViewOS,iMC) - proj2OS);
    imdiff = bprojdd(params);
    imupdateMC(:,:,iMC) = imdiff;%params.NumofView/length(indViewOS)*imdiff;
end
