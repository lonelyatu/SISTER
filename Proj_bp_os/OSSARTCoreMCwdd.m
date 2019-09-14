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

imupdateMC = iminitMC; %  ���ڵĵ������ͼ��
projMC = params.projMC; % ͶӰ����
nMC = size(projMC,3);   % �������ڵ��ؽ�ͼ��ͶӰ���ݺͱ�׼ͶӰ����֮��Ĳ��죬Ȼ��Ȼ��ͶӰ���� A_T(Ax_s - y_s)

iViews = 0 : 360/params.NumofView : 360 - 360/params.NumofView;
for iMC = 1:nMC
    
    params.iViews = iViews(indViewOS);
    params.im = imupdateMC(:,:,iMC);
    proj2OS = projdd(params);
    
    params.proj = projTotalWeightMC(:,indViewOS,iMC).*(projMC(:,indViewOS,iMC) - proj2OS);
    imdiff = bprojdd(params);
    imupdateMC(:,:,iMC) = imdiff;%params.NumofView/length(indViewOS)*imdiff;
end
