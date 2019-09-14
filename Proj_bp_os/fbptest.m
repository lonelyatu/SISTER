% ʵ�������ؽ�������
%
% 2014-06-01


addpath('D:\Yanbo\code\Spectral CT real data processing\')
reconsize = 256;
NumofView = 720;%720;%360;
NumofRou = 512;%1024;
dectsize = 0.008;%  ����õ���һ��detector bin�ĳߴ�Ϊ0.0553mm������ֵΪ0.0550mm
% FOV = 34.89;     % ɨ�跶Χֱ��
% reconreal = 18.41;
Dsourc2obj = 5;
Dsourc2detec = 10; 
Dobj2detec =  Dsourc2detec - Dsourc2obj;
pixreal = 2/256;

downsamp = 1;           % ��ʾͶӰ����2���²���
Mirror = zeros(reconsize,reconsize);
NumofTheta = NumofView;%660;
NumofRoudown = round(NumofRou/downsamp);  
dudown = dectsize*downsamp;
DofSourceCentr = Dsourc2obj;
DofSourceDetector = Dsourc2detec;

pixelsize = pixreal;
dectectorshift = 1;                     % ��ֵΪ1��ʾ�������Ԫ����˳��ߵ�������Ϊ0

YMin = -1;                              % ��4ֵΪ�ؽ�ͼ��Χ
YMax = 1;
XMin =  -1;
XMax = 1;

MagFactor = 1;

nMC = 8;
Recon = zeros(reconsize,reconsize,nMC);
for i = 1:nMC
Proj1 = reshape(flipud(Proj_1e5(:,:,i)),NumofRoudown,1,NumofView);
Recon(:,:,i) = Recon_offcenter_fullscan_pixsize_MTF_rot(Proj1,YMin,YMax,XMin,XMax,DofSourceCentr,Dobj2detec,dudown,dectsize,pixelsize,reconsize,0,360,dectectorshift);
figure,imshow(Recon(:,:,i),[]), title(num2str(i))
end
% Proj1 = reshape(flipud(fliplr(projsim(18:end-18,:))),NumofRoudown,1,NumofView);
% Reconsim = Recon_offcenter_fullscan_pixsize_MTF_rot(Proj1,YMin,YMax,XMin,XMax,DofSourceCentr,Dobj2detec,dudown,dectsize,pixelsize,reconsize,0,360,dectectorshift);
% figure,imshow(Reconsim,[])
% 
% %%
% 
% for i = 1:512
%     projbinning(i,:) = sum(proj(10*(i-1)+1:10*i, :), 1);
%         
% end




