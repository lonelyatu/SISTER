function img = imratioMC(im, imratio)
% 
% input variables:
%           'im'                  multichannel image
%           'imratio'          ratio of multichannel image
%
% Yanbo Zhang
% University of Massachusetts Lowell
% yanbozhang007@gmail.com
% 2015-07-21

nChannel = size(im,3);
img = zeros(size(im));
for i = 1:nChannel
    img(:,:,i) = im(:,:,i)*imratio(i); 
end
