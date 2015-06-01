function output = MGPtracker(N,dr,dt,pixelSize,ForwardBackward,PlotMIP)

% parameters
%   N   = number of radial segments
%   dr  = spatial step [um]
dr          = floor(dr/pixelSize);  % now converted to pixels
thresh      = 10;                   % threshold to detect the front
screenSize  = [2560 1440];          % screen pixel resolution
movieRes    = 1024;                 % output movies resolution
pos         = [(screenSize(1)-movieRes)/2 (screenSize(2)-movieRes)/2 movieRes movieRes];
poswide     = [(screenSize(1)-2*movieRes)/2 (screenSize(2)-movieRes)/2 2*movieRes movieRes];

% store input parameters
output.input.ForwardBackward    = ForwardBackward;
output.input.N                  = N;
output.input.dr                 = dr;
output.input.dt                 = dt;
output.input.pixelSize          = pixelSize;

% prompt user to select the file to analyse (for each set, select the file containting "Thresh")
[FileName,PathName,FilterIndex] = uigetfile('*.tif');

% get movie infos
nfo                         = imfinfo([PathName FileName]);
MoviesData.info             = nfo;
MoviesData.frames           = length(nfo);
MoviesData.scl              = 2^MoviesData.info(1).BitDepth-1;
MoviesData.totPixels        = MoviesData.info(1).Width*MoviesData.info(1).Height*MoviesData.frames;
clear nfo;

% get DIC image infos
nfo                         = imfinfo([PathName strrep(FileName,'Thresh','DIC')]);
DIC.info                    = nfo;
DIC.frames                  = length(nfo);
DIC.scl                     = 2^DIC.info(1).BitDepth-1;
DIC.totPixels               = DIC.info(1).Width*DIC.info(1).Height*DIC.frames;
clear nfo;

% get MIP image infos
if strcmp(PlotMIP,'yes')
    nfo                     = imfinfo([PathName strrep(FileName,'Thresh','MIP')]);
    MIP.info                = nfo;
    MIP.frames              = length(nfo);
    MIP.scl                 = 2^MIP.info(1).BitDepth-1;
    DIC.totPixels           = MIP.info(1).Width*MIP.info(1).Height*MIP.frames;
    clear nfo;
end

% initialize container for movie
if MoviesData.info(1).BitDepth == 16
    C = zeros(MoviesData.info(1).Width,...
        MoviesData.info(1).Height,...
        MoviesData.frames,...
        'uint16');
elseif MoviesData.info(1).BitDepth == 8
    C = zeros(MoviesData.info(1).Width,...
        MoviesData.info(1).Height,...
        MoviesData.frames,...
        'uint8');
else
    error('unknown bit depth');
end

% initialize container for DIC image
if DIC.info(1).BitDepth == 16
    D = zeros(DIC.info(1).Width,...
        DIC.info(1).Height,...
        DIC.frames,...
        'uint16');
elseif DIC.info(1).BitDepth == 8
    D = zeros(DIC.info(1).Width,...
        DIC.info(1).Height,...
        DIC.frames,...
        'uint8');
else
    error('unknown bit depth');
end

% initialize container for MIP image
if strcmp(PlotMIP,'yes')
    if MIP.info(1).BitDepth == 16
        mip = zeros(MIP.info(1).Width,...
            MIP.info(1).Height,...
            MIP.frames,...
            'uint16');
    elseif MIP.info(1).BitDepth == 8
        mip = zeros(MIP.info(1).Width,...
            MIP.info(1).Height,...
            MIP.frames,...
            'uint8');
    else
        error('unknown bit depth');
    end
end

% import the registered sequence
for frm = 1:1:MoviesData.frames
    C(:,:,frm) = imread([PathName FileName],'Index',frm);
end
C = C > 0;

% import the DIC sequence
D(:,:,1) = imread([PathName strrep(FileName,'Thresh','DIC')],'Index',1);

% import the MIP sequence
if strcmp(PlotMIP,'yes')
    for frm = 1:1:MIP.frames
        mip(:,:,frm) = imread([PathName strrep(FileName,'Thresh','MIP')],'Index',frm);
    end
end

% display image and ask user input for the point source of the
% chemoattractant
imshow(imadjust(D(:,:,1)));
set(gcf,'Position',pos);
title('select the point source of the chemoattractant');
[x,y] = ginput(1);
fclose all, close all;

% create patches
idx = 1;                                                            % generic index to count the patches
for i = 1:1:N
    for j = 1:dr:MoviesData.info(1).Width/2
        ang.a     = 2*pi/N*(i-1);
        ang.b     = 2*pi/N*i;
        vert{idx} = [...
            j*cos(ang.a)        j*sin(ang.a);...                    % x and y vertex coordinates
            (j+dr)*cos(ang.a)   (j+dr)*sin(ang.a);...               %
            (j+dr)*cos(ang.b)   (j+dr)*sin(ang.b);...               %
            j*cos(ang.b)        j*sin(ang.b)];                      %
        vert{idx}(:,1)  = vert{idx}(:,1)+x;                         %
        vert{idx}(:,2)	= vert{idx}(:,2)+y;                         %
        idx             = idx+1;                                    % update index
    end
end
NCircles    = length(vert)/N;
fac         = [1 2 3 4];                                            % vertices to connect

% analyse chemotaxis
if strcmp(PlotMIP,'yes')
    hdlp.image      = imshow(mip(:,:,1)); hold on;
else
    hdlp.image      = imshow(C(:,:,1)); hold on;
end
set(gcf,'Position',pos);
hdlp.crosshair  = plot(x,y,'+r');
hdlp.polygon    = [];
if strcmp(PlotMIP,'yes')
else
    for k = 1:1:length(vert)
        hdl(k) = patch('Faces',fac,'Vertices',vert{k},'FaceColor','none','EdgeColor','g');
    end
end
if exist([PathName strrep(FileName,'Thresh','MovieAnalysis')]) == 2
    delete([PathName strrep(FileName,'Thresh','MovieAnalysis')]);
end
if exist([PathName strrep(FileName,'Thresh','MIPAnalysis')]) == 2
    delete([PathName strrep(FileName,'Thresh','MIPAnalysis')]);
end
for frm = 1:1:MoviesData.frames
    delete(hdlp.image);
    if strcmp(PlotMIP,'yes')
        hdlp.image      = imshow(imadjust(mip(:,:,frm)));
    else
        hdlp.image      = imshow(C(:,:,frm));
    end
    setlayer(hdlp.image,Inf);
    for i = 1:1:N
        idx = (i-1)*NCircles+1;
        for j = 1:dr:MoviesData.info(1).Width/2
            ang.a     = 2*pi/N*(i-1);
            ang.b     = 2*pi/N*i;
            mask      = poly2mask(vert{idx}(:,1),vert{idx}(:,2),MoviesData.info(1).Width,MoviesData.info(1).Height);
            
            % avoid front going backwards
            if strcmp(ForwardBackward,'enabled')
                if frm > 1
                    if j == rgstr(i,frm-1)
                        rgstr(i,frm)        = j;
                        vertices(i,1,frm)   = (j+dr/2)*cos((ang.a+ang.b)/2)+x;
                        vertices(i,2,frm)   = (j+dr/2)*sin((ang.a+ang.b)/2)+y;
                        break
                    end
                end
            end
            
            % stop at the frame
            if      (j+dr/2)*cos((ang.a+ang.b)/2)+x > MoviesData.info(1).Width ||...
                    (j+dr/2)*sin((ang.a+ang.b)/2)+y > MoviesData.info(1).Height ||...
                    (j+dr/2)*cos((ang.a+ang.b)/2)+x < 1 ||...
                    (j+dr/2)*sin((ang.a+ang.b)/2)+y < 1
                rgstr(i,frm)        = j;
                vertices(i,1,frm)   = (j+dr/2)*cos((ang.a+ang.b)/2)+x;
                vertices(i,2,frm)   = (j+dr/2)*sin((ang.a+ang.b)/2)+y;
                break
            end
            
            % stop at the last circle
            if j >= MoviesData.info(1).Width/2-dr
                rgstr(i,frm)        = j;
                vertices(i,1,frm)   = (j+dr/2)*cos((ang.a+ang.b)/2)+x;
                vertices(i,2,frm)   = (j+dr/2)*sin((ang.a+ang.b)/2)+y;
                break
            end
            
            % detect front
            if sum(sum(C(:,:,frm).*mask)) > thresh || sum(sum(C(:,:,frm).*mask)) > sum(sum(mask))*2/3
                rgstr(i,frm)        = j;
                vertices(i,1,frm)   = (j+dr/2)*cos((ang.a+ang.b)/2)+x;
                vertices(i,2,frm)   = (j+dr/2)*sin((ang.a+ang.b)/2)+y;
                break
            end
            
            % update index
            idx     = idx+1;
        end
        rng(i,frm)  = sqrt((vertices(i,1,frm)-x)^2+(vertices(i,2,frm)-y)^2);
    end
    delete(hdlp.polygon);
    hdlp.polygon    = plot([vertices(:,1,frm);vertices(1,1,frm)],[vertices(:,2,frm);vertices(1,2,frm)],'-r','LineWidth',2);
    area(frm)       = sum(sum(poly2mask(vertices(:,1,frm),vertices(:,2,frm),MoviesData.info(1).Width,MoviesData.info(1).Height)));
    dist(frm)       = sum(rng(:,frm))/N;
    pause(0.01)
    if strcmp(PlotMIP,'yes')
        imwrite(frame2im(getframe(gcf,[0.11 0.11 0.78 0.78]*movieRes)),[PathName strrep(FileName,'Thresh','MIPAnalysis')],'WriteMode','append');
    else
        imwrite(frame2im(getframe(gcf,[0.11 0.11 0.78 0.78]*movieRes)),[PathName strrep(FileName,'Thresh','MovieAnalysis')],'WriteMode','append');
    end
end

% create output structure
output.area         = area*pixelSize^2; % converts to um^2
output.area_norm	= area/max(area);   % normalised area
output.x            = x;
output.y            = y;
output.dist         = dist*pixelSize;   % converts to um

% calculate 5% cut-off time
[M,I] = max(output.area_norm < 0.05);
output.area_cutoff = (I-1)*dt;
[M,I] = max(output.dist/max(output.dist) < 0.05);
output.dist_cutoff = (I-1)*dt;

if frm > 1
    % Fit line through area
    coeffNames      = {'a','b','c'};
    myfun   = fittype(...
        'a*exp(-x/b)+c',...
        'independent','x',...
        'coefficients',coeffNames);
    options = fitoptions(...
        'method','NonLinearLeastSquares',...
        'StartPoint',[100 5 0],...
        'MaxFunEvals',10000,...
        'TolFun',1e-07,...
        'TolX',1e-07,...
        'Lower',[0 0 0],...
        'Upper',[+Inf +Inf +Inf]);
    timeScale = dt*(0:1:MoviesData.frames-1);
    [cfun_area,gof_area] = fit(timeScale',output.area',myfun,options);
    
    % Fit line through mean distance
    [cfun_dist,gof_dist] = fit(timeScale',output.dist',myfun,options);
    
    % figure for area and mean distance
    fig = figure;
    set(fig,'Position',poswide);
    subplot(121);
    plot(timeScale,output.area,'-k'); hold on;
    plot(timeScale,cfun_area(timeScale),'-r','LineWidth',2);
    text(max(timeScale)/2,1.1*max(output.area)/2,['\tau = ' num2str(cfun_area.b,'%4.1f') ' min']);
    xlabel('Time (min)');
    ylabel('Clear area (µm^2)');
    subplot(122);
    plot(timeScale,output.dist,'-k'); hold on;
    plot(timeScale,cfun_dist(timeScale),'-r','LineWidth',2);
    text(max(timeScale)/2,1.1*max(output.dist)/2,['\tau = ' num2str(cfun_dist.b,'%4.1f') ' min']);
    xlabel('Time (min)');
    ylabel('Average distancearea (µm)');
    
    % update output structure
    output.fit_area     = cfun_area;
    output.fit_dist     = cfun_dist;
    
    % save figure as a JPEG
    set(fig,...
        'PaperUnits','inches',...
        'PaperSize',[8 4],...
        'PaperPosition',[0 0 8 4]);
    print(fig,[PathName strrep(strrep(FileName,'Thresh','Output'),'.tif','.jpeg')],'-djpeg','-r300');
end

% save output
if exist([PathName strrep(strrep(FileName,'Thresh','Output'),'.tif','.m')]) == 2
    delete([PathName strrep(strrep(FileName,'Thresh','Output'),'.tif','.m')]);
end
save([PathName strrep(strrep(FileName,'Thresh','Output'),'.tif','.mat')],'output');

% the end