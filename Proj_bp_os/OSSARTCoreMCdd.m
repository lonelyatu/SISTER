% 基于distance driven前后向投影的OS-SART核心
% 针对多通道数据
%
% 2015-10-14

function imupdateMC = OSSARTCoreMCdd(params, iminitMC, ios, OSnum)

% iternum = 20;
% OSnum = 20;
imupdateMC = iminitMC;
projMC = params.projMC;
nMC = size(projMC,3);

iViews = 0 : 360/params.NumofView : 360 - 360/params.NumofView;
% for ios = 1:OSnum
for iMC = 1:nMC
    indViewOS = ios:OSnum:params.NumofView;
    
    params.iViews = iViews(indViewOS);
    params.im = imupdateMC(:,:,iMC);
    proj2OS = projdd(params);
    
    params.proj = projMC(:,indViewOS,iMC) - proj2OS;
    imdiff = bprojdd(params);
    imupdateMC(:,:,iMC) = imdiff;%params.NumofView/length(indViewOS)*imdiff;
end
% end
