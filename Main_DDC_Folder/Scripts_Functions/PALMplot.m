function  Density_Image=PALMplot(coords, coords2, suffix)
    % 'coords' has three columns: [x y sigma] in nm
    % if coords is only two columns, the software will use the default
    % 'sigma_precision' value instead
    % 'suffix' is an optional string to add to the filename
    % if your coordinates are in pixels, change 'Pixel_sz' to the pixel size in
    % nm (sigma still assumed to be in in nm)
    
    % edited 5/19/2015 CC to fix some small rounding errors in defining the
    % center of each molecule's plotted Gaussian.

    % plotting parameters
    Pixel_sz = 1;           % nm, assume coords in nm
    sigma_precision = 40 ;  % nm, default localization precision if none provided    
    sigmaThresh = Inf;      % use this to filter out spots with large sigmas
    PALMPixelSize = 10;     % nm, output PALM image pixel size   
    
    % add localization precision column if required
    if size(coords,2) < 3
        coords(:,3) = sigma_precision;
        coords2(:,3) = sigma_precision;
    end
    
    % name suffix
    if nargin < 2
        suffix = '';
    end
        
    % Gaussian calucation params
    GaussPlotRange = 21;      % number of PALM pixels wide to calculate each Gaussian spot
    r = (GaussPlotRange - 1) / 2; 
    pixMultiple  = Pixel_sz / PALMPixelSize;

    % generate meshgrid for simulated 2D gaussians
    GaussPlotCenter = (GaussPlotRange + 1) / 2;
    [xPALM, yPALM] = meshgrid(1:GaussPlotRange);

    
%    coords((coords(:,2)>.5*10^4),:)=[];
%    coords((coords(:,1)>.5*10^4),:)=[];
    
    % calculate necessary image size        
    
    
    %CHB, I am modifying this so that all of the figures are the same size
    frl = floor(min(coords2(:,2)));  % row, low
    frh = ceil(max(coords2(:,2)));   % row, high
    fcl = floor(min(coords2(:,1)));  % column, low
    fch = ceil(max(coords2(:,1)));   % column, high
    nRows = frh - frl + 1;  
    nCols = fch - fcl + 1;  
        
    % rescaled PALM image size
    zoomRows = round(nRows * pixMultiple) + 4*r; % 2r border on each edge
    zoomCols = round(nCols * pixMultiple) + 4*r; % 2r border on each edge
           
    % extract coordinate/sigma data
    devX  = coords(:,3);
    devY  = coords(:,3);
    xx0c  = coords(:,1) * pixMultiple;
    yy0c  = coords(:,2) * pixMultiple;
    
    xx0c2  = (coords2(:,1)) * pixMultiple;
    yy0c2  = (coords2(:,2)) * pixMultiple;
    
    % zero x and y values
    xx0c = xx0c - min(xx0c2);
    yy0c = yy0c - min(yy0c2);
         
    % allocate image
    palmImg = double(zeros(zoomRows,zoomCols));
    %size( palmImg)
    % loop through each molecule
    nMolec = length(xx0c);
    for iMolec = 1:nMolec
        size( palmImg);
       
        % if precision is good enough, add to palm image
        if devX(iMolec) <= sigmaThresh && devY(iMolec) <= sigmaThresh
        
            % molecule position
            yPos = round(yy0c(iMolec)) + 2*r; % 2r border on each edge
            xPos = round(xx0c(iMolec)) + 2*r; % 2r border on each edge   
            deta2x = devX(iMolec)/ PALMPixelSize;   % loc. precision in PALM pixels
            deta2y = devY(iMolec)/ PALMPixelSize;   % loc. precision in PALM pixels

            % if molecule is within image bounds
            if ~(yPos-r<1 || yPos+r>zoomRows || xPos-r<1 || xPos+r>zoomCols)
                 
                % create gaussian spot
                yCenter = yy0c(iMolec) + 2*r - yPos;  % distance of real center from center pixel
                xCenter = xx0c(iMolec) + 2*r - xPos;  % distance of real center from center pixel   
                expNumeratorX = -(xPALM-(GaussPlotCenter+xCenter)).^2;
                expNumeratorY = -(yPALM-(GaussPlotCenter+yCenter)).^2;           
                spot_gaussian = 1/(2*pi*deta2x*deta2y)*exp(expNumeratorX/(2*deta2x.^2) + expNumeratorY/(2*deta2y.^2));
                
                % add gaussian spot to palm images
                palmImg(yPos-r:yPos+r, xPos-r:xPos+r) = ...
                palmImg(yPos-r:yPos+r, xPos-r:xPos+r) + spot_gaussian(1:2*r+1,1:2*r+1);

            end  % out of bounds condition     

        end
        
    end  % end iMolec loop    
          
    % write palm image
    %{
    if max(palmImg(:)) < 256
        imwrite(uint16(palmImg/max(palmImg(:))*256),['palm_' suffix '.tif'],'tif','Compression','none','WriteMode','overwrite');
    else
        imwrite(uint16(palmImg),['palm_' suffix '.tif'],'tif','Compression','none','WriteMode','overwrite');
    end
    %}
    % print mean sigma
   % disp(['Mean sigma: ' num2str(mean([devX; devY])) ' nm'])
    
    % print number of molecules
   % disp(['# of molecules: ' num2str(nMolec)])
    
    % show palm image
    %figure(sum('PALMplotN'))
    %imagesc(palmImg); axis equal; axis off; colormap(hot); set(gcf,'Name',pwd,'NumberTitle','off');
    
    
    Density_Image=palmImg;
    
    
end
