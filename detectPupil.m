function [B,R,pupilCenter,Res,imDiv] = detectPupil(img,imDiv,pupilCenterPrev,pupilBounds,nThBins,stepSize)

% startburst like detection algorithm


%%
%     img = ImGray;
%     imDiv = [];
%     pupilCenterPrev = [240 320];
%     pupilBounds = [50 160];
%     pupilCenterPrev = [540 960];
%     pupilBounds = [80 350];
%     nThBins = 16;
%     stepSize = 5;

    stepMultFact = 0.5;
    
    showPlots = 0;
    
    pixRes = 1/4; % resolution of radius annuli; 1/X; the higher X, the faster

    imSize = size(img);
    img = double(img);
    
    nRadBins = floor((pupilBounds(2)-pupilBounds(1))*pixRes);
    thsBinEdges = linspace(-pi,pi,nThBins+1);
    radsBinEdges = linspace(pupilBounds(1),pupilBounds(2),nRadBins+1);
    
    if length(imDiv) == 0
        
        X = repmat(1:imSize(2),imSize(1),1)-round(pupilCenterPrev(2));
        Y = repmat([1:imSize(1)]',1,imSize(2))-round(pupilCenterPrev(1));
        [TH,RAD] = cart2pol(X,Y);
        imDiv = zeros(imSize);
        
        countDiv = 0;
        for i = 1:nThBins
            thVect = TH >= thsBinEdges(i) & TH < thsBinEdges(i+1);
    %             subplot(6,6,i)
    %             imagesc(thVect)
    
            for j = 1:nRadBins
                countDiv = countDiv+1;
                radVect2 = RAD>radsBinEdges(j)&RAD<radsBinEdges(j+1);
                
                imDiv(radVect2&thVect) = countDiv;
            end

        end
    end
    if showPlots
    figure();
    subplot(2,2,1);
    imshow(img,[]);
    subplot(2,2,2);
    imshow(imDiv,[]); % image with wedge divisions
    end
    
    measurements = regionprops(imDiv,img,'MeanIntensity');
    if length(measurements) == nRadBins*nThBins
        imLumPerTH_RAD = reshape([measurements.MeanIntensity],[nRadBins nThBins])';
    else  % this happens if pupil is detected near the border ...
        % check whether parts in the iris are missing (because pupil is too
        % close to the image edge)
        imDivIdx = unique(imDiv);
        nImDivs = sum(imDivIdx>0);
        imDivVect = zeros(1,nRadBins*nThBins);
        imDivVect(imDivIdx(imDivIdx>0)) = 1;
        imDivVect = logical(imDivVect);
    
        tempMeanIntensity = NaN(1,nRadBins*nThBins);
        tempMeanIntensity2 = [measurements.MeanIntensity];
        
        tempMeanIntensity(imDivVect) = tempMeanIntensity2(isfinite(tempMeanIntensity2));
        imLumPerTH_RAD = reshape(tempMeanIntensity,[nRadBins nThBins])';
    end
    
    
%% increase resolution 

    imLumPerTH_RAD_hr = NaN([size(imLumPerTH_RAD,1) size(imLumPerTH_RAD,2)/pixRes]);
    imLumPerTH_RAD_hr(:,1:(1/pixRes):end) = imLumPerTH_RAD;
    imLumPerTH_RAD_hr(:,end) = imLumPerTH_RAD(:,end);
    
    imLumPerTH_RAD_hr = imLumPerTH_RAD_hr';
    finIdx = find(isfinite(imLumPerTH_RAD_hr(:)));
    intIdx = find(~isfinite(imLumPerTH_RAD_hr(:)));
    imLumPerTH_RAD_hr(intIdx) = interp1(finIdx,imLumPerTH_RAD_hr(finIdx),intIdx);
    imLumPerTH_RAD_hr = imLumPerTH_RAD_hr';
    
    if showPlots
    figure();
%     imshow(imLumPerTH_RAD_hr,[]); % rows = angles, cols = radii
    imagesc(imLumPerTH_RAD_hr); % rows = angles, cols = radii
    ylabel('Angles');
    xlabel('Radii');
    a = colorbar;
    a.Label.String = 'Luminance RGB';
    title('Higher resolution');
    end
    
%%  detect edge with luminance change threshold

%     xdata = radsBinEdges(1:stepSize:end-1);
%     xdata = xdata(1:end-1)+diff(xdata);
    xdata = linspace(pupilBounds(1),pupilBounds(2),(nRadBins/pixRes)-stepSize);
    
    
%     allDiff     = NaN(nThBins,length(xdata));
%     B           = NaN(nThBins,2);
%     RperTH      = NaN(nThBins,1);
%     
%     for i = 1:nThBins
%         tempRadArr = imLumPerTH_RAD_hr(i,:);
%         allDiff(i,:) = diff([tempRadArr(1:end-stepSize); tempRadArr(stepSize+1:end)]);
% %         allDiff(i,:) = diff(imLumPerTH_RAD_hr(i,1:stepSize:end)')';
%         
%         rIdx = find(allDiff(i,:)>stepSize,1,'first');
%         if rIdx
%             [pks,locs] = findpeaks(allDiff(i,rIdx:end));
%             if length(locs) == 0 % maximum
%                 [~,locs] = max(allDiff(i,rIdx:end));
%             end
%             RperTH(i) = xdata(locs(1)+rIdx-1);
%             [B(i,2),B(i,1)] = pol2cart(mean([thsBinEdges(i) thsBinEdges(i+1)]),RperTH(i));
%         end
%     end
    
    allDiff     = NaN(nThBins,length(xdata));
    B           = NaN(nThBins,2);
    RperTH      = NaN(nThBins,1);
    
    for i = 1:nThBins
        tempRadArr = imLumPerTH_RAD_hr(i,:);
        allDiff(i,:) = diff([tempRadArr(1:end-stepSize); tempRadArr(stepSize+1:end)]);
%         allDiff(i,:) = diff(imLumPerTH_RAD_hr(i,1:stepSize:end)')';
        
        rIdx = find(allDiff(i,:)>stepSize*stepMultFact,1,'first')-1;
        if rIdx
            [pks,locs] = findpeaks(allDiff(i,rIdx:end));
            if length(locs) == 0 % no peaks detected
                [~,locs] = max(allDiff(i,rIdx:end)); % take maximum
                radIdx = locs(1)+rIdx-1;
                RperTH(i) = xdata(radIdx);
%             elseif length(locs) > 1 % more than 1 peak detected
%                 [~,locsIdx] = max(allDiff(i,locs+rIdx-1));
%                 radIdx = locs(locsIdx(1))+rIdx-1;
%                 RperTH(i) = xdata(radIdx);
            else
                radIdx = locs(1)+rIdx-1;
                RperTH(i) = xdata(radIdx);
            end
        else
            [~,locs] = max(allDiff(i,:));
            radIdx = locs;
            RperTH(i) = xdata(radIdx);
        end
        [B(i,2),B(i,1)] = pol2cart(mean([thsBinEdges(i) thsBinEdges(i+1)]),RperTH(i));
        
    end

    
    
    %% x,y coordinates of polygon pupil border 

    B = B(isfinite(B(:,1)),:);    
    B(:,1) = B(:,1)+pupilCenterPrev(1);
    B(:,2) = B(:,2)+pupilCenterPrev(2);
    
    %% extract properties
    
    imBW = poly2mask(B(:,2),B(:,1),imSize(1),imSize(2));
    imBW = bwconvhull(imBW,'objects'); % fill in concave shape of pupil (sometimes pupil is not fully filled)
    
%     Also make convex
    measurements = regionprops(imBW,'Circularity','Area','Centroid','MajorAxisLength','MinorAxisLength');
    
    if sum(imBW(:)) > pi*(pupilBounds(1)^2)
        if length(measurements) > 1 % select the one with the largest area
            allAreas = [];
            for i = 1:length(measurements)
                allAreas(i) = measurements(i).Area;
            end
            [~,maxIdx] = max(allAreas);
            measurements = measurements(maxIdx);
        end
        
        Res = abs(measurements.Circularity-1);
        centers = measurements.Centroid;
        diameters = mean([measurements.MajorAxisLength measurements.MinorAxisLength],2);
        R = diameters/2;
        pupilCenter = fliplr(centers);
        
    else
        Res = 999;
        R = NaN;
        pupilCenter = pupilCenterPrev;
        B = [NaN NaN];
    end


%% plot results

    if showPlots
        figure();
        plot(radsBinEdges(1:end-1),imLumPerTH_RAD');
        xlabel('Radii');
        ylabel('Luminance');

        colormapper = colormap('jet');
        colormapper = colormapper(round(linspace(1,256,size(allDiff,1))),:);

        figure(99);
        for i = 1:size(allDiff,1)
            plot(xdata,allDiff(i,:),'k','Color',colormapper(i,:));
            hold on
            rIdx = find(allDiff(i,:)>stepSize*stepMultFact,1,'first')-1;
            if rIdx
                [pks,locs] = findpeaks(allDiff(i,rIdx:end));
                if length(locs) == 0 % maximum
                    [~,locs] = max(allDiff(i,rIdx:end));
                    radIdx = locs(1)+rIdx-1;
                    RperTH(i) = xdata(radIdx);
%                 elseif length(locs) > 1 % more than 1 peak detected
%                     [~,locsIdx] = max(allDiff(i,locs+rIdx-1));
%                     radIdx = locs(locsIdx(1))+rIdx-1;
%                     RperTH(i) = xdata(radIdx);
                else
                    radIdx = locs(1)+rIdx-1;
                    RperTH(i) = xdata(radIdx);
                end
            else
                [~,locs] = max(allDiff(i,:));
                radIdx = locs;
            end

            plot(xdata(radIdx),allDiff(i,radIdx),'ko','Color',colormapper(i,:),'MarkerSize',20,'MarkerFaceColor',colormapper(i,:))
            text(xdata(radIdx),allDiff(i,radIdx),num2str(i),'Color','w','HorizontalAlignment','center')

        end
        plot(xdata,nanmean(allDiff),'k','LineWidth',4);
        line([xdata(1) xdata(end)],[stepSize*stepMultFact stepSize*stepMultFact],'Color','k','LineStyle',':','LineWidth',4)
        hold off
        xlabel('Radii');
        ylabel('Diff. in Luminance across Radii');

        figure(999);
        imshow(uint8(img));
        hold on
        [xx,yy] = pol2cart(linspace(-pi,pi,nThBins),repmat(R,1,nThBins));
        
        plot(pupilCenter(2),pupilCenter(1),'go');
        for i = 1:length(xx)
            plot(xx(i)+pupilCenter(2),yy(i)+pupilCenter(1),'g*');
            text(xx(i)+pupilCenter(2),yy(i)+pupilCenter(1),num2str(i),'Color','w')
        end
        
        plot(pupilCenterPrev(2),pupilCenterPrev(1),'bo');
        plot(xx+pupilCenterPrev(2),yy+pupilCenterPrev(1),'bx');

        plot(mean(B(:,2)),mean(B(:,1)),'ro');
        plot(B(:,2),B(:,1),'r+');

        title('Red = fitted area; Blue = previous border; Green = new border')

    end
    
