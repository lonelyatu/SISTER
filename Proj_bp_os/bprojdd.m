function im = bprojdd(params)
% Distance driven backprojection for equi-distance fanbeam CT
% Reference:
% B. De Man and S. Basu, "Distance-driven projection and
% backprojection in three dimensions," Physics in Medicine
% and Biology, vol. 49, pp. 2463-2475, 2004.
%
% input field of params:
%           'proj'                      projection
%           'Dsource2centr'     distance from source to center, mm
%           'Dsource2detec'    distance from source to detector, mm
%           'pixelsize'              pixel size, mm
%           'binsize'                 detector bin size, mm
%           'binshift'                detector bin shift, mm
%           'iViews'                  view angels, degree, optional
%
%

%% fanbeam parameters

proj = params.proj;
reconsize = params.reconsize;
Dsource2centr = params.Dsource2centr;
Dsource2detec = params.Dsource2detec;
[NumofBin, NumofView] = size(proj);
pixelsize = params.pixelsize;
binsize = params.binsize;
binshift = params.binshift;               % detector shift, mm

% CT scanning views
if isfield(params, 'iViews')
    iViews = params.iViews;
else
    startAngle = 0;
    endAngle = 360-360/NumofView;
    iViews = startAngle : (endAngle - startAngle)/(NumofView - 1) : endAngle;               % degree
end

% if CT scanning is clockwise
if isfield(params, 'clockwise')
    if params.clockwise == 1
        iViews = mod(360 - iViews, 360);
    end
end

iBinPos = binsize*(-NumofBin/2:NumofBin/2) + binshift;    % rays' boudaries, NumofBin+1, mm
im = zeros(reconsize);

%% main, equi-distance detector

for ip =  1:length(iViews)
    iview = iViews(ip);
    ibeta = mod(iview/180*pi, 2*pi);
    
    % position of x-ray source
    xSource = -Dsource2centr*sin(ibeta);
    ySource = Dsource2centr*cos(ibeta);
    
    % projection for current view
    projline = proj(:, ip);
    
    % source near y axis, so compute along x axis
    if (ibeta < pi/4) || (abs(ibeta - pi) < pi/4) || (2*pi - ibeta< pi/4)
        % map all the rays' boundaries to x axis
        gammas = atan(iBinPos/Dsource2detec) ;
        detecPos = xSource + ySource*tan(gammas + ibeta);
        if  abs(ibeta - pi) < pi/4
            % in this case, the position of detector is decreasing, we need
            % to reverse it
            detecPos = detecPos(end:-1:1);
            projline = projline(end:-1:1);
        end
        deteclen = detecPos(2:end) - detecPos(1:end-1);
        % cos of projection line
        cosRay = ySource./sqrt((detecPos - xSource).^2 + ySource^2);
        cosRayMid = abs((cosRay(1:end-1)+cosRay(2:end))/2)'; % abs ray cos
        projlineWeight = projline./cosRayMid;       % weighted projection
        
        % map all the image pixel boundaries to x axis
        linebound = ([0:reconsize] - reconsize/2);
        imbound = repmat(linebound, reconsize, 1);
        imy = repmat([1:reconsize]' - reconsize/2 - 0.5, 1, reconsize+1); % y pos for each line pixels
        imPos = xSource + ySource*(imbound*pixelsize - xSource)./(ySource - imy*pixelsize);
        
        for ypix = 1:reconsize
            
            imline = zeros(1, reconsize);
            ix = 2;
            id = 2;
            xPosLine = imPos(ypix,:);
            currentPos = max([detecPos(1), xPosLine(1)]); % find the start position
            dPos = detecPos(id);
            xPos = xPosLine(ix);
            
            % if xPos <= currentPos, find start postion of image
            while xPos <= currentPos
                ix = ix + 1;
                if ix > reconsize+1
                    break;
                else
                    xPos = xPosLine(ix);
                end
            end
            
            % if dPos <= currentPos, find start position of detector
            while dPos <= currentPos
                id = id + 1;
                if id > NumofBin+1
                    break;
                else
                    dPos = detecPos(id);
                end
            end
            
            % main loop
            while ix <= reconsize + 1 && id <= NumofBin + 1
                xPos = xPosLine(ix);
                dPos = detecPos(id);
                id_1 = id - 1;
                ix_1 = ix - 1;
                if xPos < dPos
                    len = xPos - currentPos;
                    imline(ix_1) =  imline(ix_1) + len/deteclen(id_1)*projlineWeight(id_1);
                    currentPos = xPos;
                    ix = ix + 1;
                elseif xPos > dPos
                    len = dPos - currentPos;
                    imline(ix_1) =  imline(ix_1) + len/deteclen(id_1)*projlineWeight(id_1);
                    currentPos = dPos;
                    id = id + 1;
                else
                    len = xPos - currentPos;
                    imline(ix_1) =  imline(ix_1) + len/deteclen(id_1)*projlineWeight(id_1);
                    currentPos = xPos;
                    ix = ix + 1;
                    id = id + 1;
                end
            end
            
            im(ypix, :) = im(ypix, :) + imline/pixelsize;
        end
        
    else
        
        % map all the rays' boundaries to y axis
        gammas = atan(iBinPos/Dsource2detec) ;
        detecPos = ySource + xSource*tan(pi/2-(gammas + ibeta));
        if  abs(ibeta - 3/2*pi) <= pi/4
            % in this case, the position of detector is decreasing, we need
            % to reverse it
            detecPos = detecPos(end:-1:1);
            projline = projline(end:-1:1);
        end
        deteclen = detecPos(2:end) - detecPos(1:end-1);
        % cos of projection line
        cosRay = xSource./sqrt((detecPos - ySource).^2 + xSource^2);
        cosRayMid = abs((cosRay(1:end-1)+cosRay(2:end))/2)'; % abs ray cos
        projlineWeight = projline./cosRayMid;       % weighted projection
        
        % map all the image pixel boundaries to y axis
        linebound = ([0:reconsize] - reconsize/2)';
        imbound = repmat(linebound, 1, reconsize);
        imx = repmat([1:reconsize] - reconsize/2 - 0.5, reconsize+1, 1); % x pos for each line pixels
        imPos = ySource + xSource*(imbound*pixelsize - ySource)./(xSource - imx*pixelsize);
        
        for xpix = 1:reconsize
            
            imline = zeros(reconsize, 1);
            iy = 2;
            id = 2;
            yPosLine = imPos(:, xpix);
            currentPos = max([detecPos(1), yPosLine(1)]); % find the start position
            dPos = detecPos(id);
            yPos = yPosLine(iy);
            
            % if xPos <= currentPos, find start postion of image
            while yPos <= currentPos
                iy = iy + 1;
                if iy > reconsize+1
                    break;
                else
                    yPos = yPosLine(iy);
                end
            end
            
            % if dPos <= currentPos, find start position of detector
            while dPos <= currentPos
                id = id + 1;
                if id > NumofBin+1
                    continue;
                else
                    dPos = detecPos(id);
                end
            end
            
            % main loop
            while iy <= reconsize + 1 && id <= NumofBin + 1
                yPos = yPosLine(iy);
                dPos = detecPos(id);
                id_1 = id - 1;
                iy_1 = iy - 1;
                if yPos < dPos
                    len = yPos - currentPos;
                    imline(iy_1) =  imline(iy_1) + len/deteclen(id_1)*projlineWeight(id_1);
                    currentPos = yPos;
                    iy = iy + 1;
                elseif yPos > dPos
                    len = dPos - currentPos;
                    imline(iy_1) =  imline(iy_1) + len/deteclen(id_1)*projlineWeight(id_1);
                    currentPos = dPos;
                    id = id + 1;
                else
                    len = yPos - currentPos;
                    imline(iy_1) =  imline(iy_1) + len/deteclen(id_1)*projlineWeight(id_1);
                    currentPos = yPos;
                    iy = iy + 1;
                    id = id + 1;
                end
            end
            
            im(:, xpix) = im(:, xpix) + imline/pixelsize;
        end
        
    end
    
end