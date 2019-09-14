% OS-SARTÀ„∑®
%
% 2015-10-14

function imupdate = OSSART(params, iminit, ios, OSnum, bprojones)

% iternum = 20;
% OSnum = 20;
imupdate = iminit;
proj = params.proj;

iViews = 0 : 360/params.NumofView : 360 - 360/params.NumofView;
% for ios = 1:OSnum
    indViewOS = ios:OSnum:params.NumofView;
    
    params.iViews = iViews(indViewOS);
    params.im = imupdate;
    proj2OS = projdd(params);
    
    params.proj = proj(:,indViewOS) - proj2OS;
    imdiff = bprojdd(params);
    imupdate = imupdate + params.NumofView/length(indViewOS)*imdiff./bprojones;
% end