function [eyeData,eyeDataHeader,vidInfo] = irissometry(vidFileName,configVar)
% eyeData = irissometry(vidFileName,configVar)
%
% -------------------
% Input
% -------------------
% 
% vidFileName: string with full path and filename to a video.
% 
% configVar: struct with the following variables for irissometry 
% configuration:
% 
% .scleraRadius: integer indicating the radius of the limbus (border between
% the iris and sclera (eye-white) in pixels. 
% 
% .pupilRadiusRange: array with two integers. Integers indicate minimum and
% maximum possible pupil size to speed up and ease 
% search for pupil borders
% 
% .startFrame: integer indicating which video frame to start the analysis.
% Increase this number if your camera needs some frames to dynamically
% adjusts luminance (always inspect first couple frames to check).
% 
% .minCircleFitResidual: float indicating the maximum allowed error. 
% 
% .nFeaturePointRange:array with two integers. Integers indicate minimum and
% maximum possible features to track in iris. To few (<300), reduces 
% accuracy of feature distance calculations. Too many (>2000), reduces 
% video processing speed to unacceptable rates.
% 
% 
% .waitForNFramesAfterPupilMiss: integer indicating for how many frames to
% wait after a missing pupil (e.g., caused by a blink) before features are
% re-initialized. Sometimes the pupil is re-detected while the lid is still
% slightly over the pupil. Wait from that moment on for the eyelid to be
% fully recovered to eyes-open position, otherwise no features are tracked
% in the upper region of the iris.
%  
% overwriteOutputData: boolean to indicate to overwrite already existing
% data output (videoFileName_output.mat).
% 
% .starburst.nThBins: integer indicating number of wedges extending from 
% pupil center to iris to calculate sudden luminance increases for pupil 
% border detection. Too few (<4) causes bad pupil border fits because not
% enough border points. Too many (>24), reduces video processing speed to 
% unacceptable rates.
% 
% .startburst.stepSize: integer indicating distance in pixels for sliding
% window to compare difference in luminance values. 5 is recommended for
% 640 pixel width video frames. Multiple this number depending on the
% increase in your videos resolution in ratio to 640 width (e.g., a video
% 1280 pixels in width should have a stepsize of 10).
% 
% .disp.saveVideoIrissometryOutput: boolean (1 = True, or 0 = False) to 
% indicate whether or not to save a bew video of the irissometry output on
% top of the original video frames
% 
% .disp.inspectFrames: boolean to indicate whether or not to show
% each video frame with irissometry output
% 
% .disp.inspectFramesManualPace: boolean to indicate whether or 
% not to show each video frame but only continues to next frame after
% SPACEBAR is pressed.
% 
% 
% Example of creating a configuration file:
% configVar = struct();
% configVar.scleraRadius = 300;
% configVar.pupilRadiusRange = [50 160];
% configVar.startFrame = 1;
% configVar.maxCircleFitError = 0.16;
% configVar.nFeaturePointRange = [600 2000];
% configVar.waitForNFramesAfterPupilMiss = 10;
% configVar.overwriteOutputData = 0;
% configVar.starburst.nThBins = 16;
% configVar.startburst.stepSize = 5;
% configVar.disp.inspectFrames = 1;
% configVar.disp.inspectFramesManualPace = 1;
% configVar.disp.saveVideoIrissometryOutput = 1;
% 
% 
% -------------------
% Output
% -------------------
%
% The script automatically saves a file with the following name:
% "videoFileName_output.mat".
% It automatically reads this if already existing. If you prefer to
% overwrite existing files, set configVar.overwriteOutputData to 1 (True)
% 
% eyeData: matrix with a row for each video frame and columns for
% the following measurements:
% 
% 1 = pupil area
% 2 = pupil radius
% 3 = x horizontal position of pupil center
% 4 = y vertical position of pupil center
% 5 = error of circle fit to pupil border (circularity of pupil). If low, indicates blink
% 6 = video frame timestamps
% 7 = number of points (wedges) detected on the pupil border
% 8 = number of feature points detected in iris
% 9 = mean distance between all iris feature points
% 
% 10:19 = iris feature distance density within a selected annulus within the iris
% For each column increase, starting at 20, the borders of the annulus move 
% to more eccentric (peripheral) regions, from 0-10% to 90-100% in steps of
% 10%.
% 
% 20:29 = number of tracked SURF points per annulus
% 
% 30 = X coordinates of fitted circle to all features
% 31 = X coordinates of fitted circle to all inner iris features
% 32 = X coordinates of fitted circle to all outer iris features
% 33:42 = X coordinates of fitted circle to features per annulus
% 
% 43 = Y coordinates of fitted circle to all features
% 44 = Y coordinates of fitted circle to all inner iris features
% 45 = Y coordinates of fitted circle to all outer iris features
% 46:55 = Y coordinates of fitted circle to features per annulus

eyeDataHeader = {'Pupil area','Pupil radius','Pupil X','Pupil Y','Pupil fit error','Timestamp','# Pupil points','# Iris points (Ip)','Ip radial distance (rd)'};
for i = 1:10
    eyeDataHeader{9+i} = ['Ip rd annu. ' num2str((i-1)*10) '-' num2str(i*10)];
    eyeDataHeader{19+i} = ['# Ip annu. ' num2str((i-1)*10) '-' num2str(i*10)];
end
eyeDataHeader{30} = 'Ip X';
eyeDataHeader{31} = 'Ip inner rim X';
eyeDataHeader{32} = 'Ip outer rim X';
eyeDataHeader{43} = 'Ip Y';
eyeDataHeader{44} = 'Ip inner rim Y';
eyeDataHeader{45} = 'Ip outer rim Y';

for i = 1:10
    eyeDataHeader{32+i} = ['Ip annu. X ' num2str((i-1)*10) '-' num2str(i*10)];
    eyeDataHeader{45+i} = ['Ip annu. Y' num2str((i-1)*10) '-' num2str(i*10)];
end

if configVar.overwriteOutputData | ~exist([vidFileName '_output.mat'])

    disp(['Analyzing: ' vidFileName '_output.mat'])
    
    % Create the point tracker object.
    pointTracker = vision.PointTracker('MaxBidirectionalError', 2);

    % load video for analysis
    vidObj = VideoReader(vidFileName);

    % save video frame-by-frame with irissometry output
    if configVar.disp.saveVideoIrissometryOutput
        vidFileName_fit = [vidFileName '_output.mp4'];
        vidObj_fit = VideoWriter(vidFileName_fit,'MPEG-4');
        vidObj_fit.Quality = 95;
        vidObj_fit.FrameRate = vidObj.FrameRate;
        open(vidObj_fit)
    end

            
    % create struct to save video information in            
    vidInfo = struct();
    vidInfo.nFrames         = vidObj.NumberOfFrames;
    vidInfo.vidDuration     = vidObj.Duration;
    vidInfo.frameRate       = (vidObj.NumberOfFrames-1)/vidInfo.vidDuration;
    vidInfo.frameRate_alt   = vidObj.FrameRate;
    vidInfo.frameNums_sel   = configVar.startFrame:vidInfo.nFrames;
    vidInfo.nFrames_sel     = length(vidInfo.frameNums_sel);
    vidInfo.tStamp          = zeros(1,vidInfo.nFrames_sel);


    % Output info about the video
    disp(['Number of frames: ' num2str(vidInfo.nFrames)]);
    disp(['Frame rate: ' num2str(vidInfo.frameRate)]);
    disp(['Frame rate alternative: ' num2str(vidInfo.frameRate_alt)]);
    disp(['Video duration: ' num2str(vidInfo.vidDuration)]);
    disp(['-------------------------------------------------------']);
    disp(['-------------------------------------------------------']);

    numPts                 = 0;        % number of points detected
    points                 = [];       % x and y coordinates of points
    nTimesPointsRedetected = 0;

    pupilDetected          = 0;
    pointsDetected         = 0;

    nFramesSeqPupilDetected  = 0;
    nFramesNoPupilDetected   = 0;



    tic
    countF = 0;
    for f = configVar.startFrame:round(vidInfo.nFrames_sel) % loop through each image frame to calculate average pixel value in region of interest / 1400:round(vidInfo.nFrames_sel)%
        countF  = countF +1;
        if sum(f == round(linspace(1,round(vidInfo.nFrames_sel),11)))
            disp(['Extracting pixel values in video: ' num2str(round((f/vidInfo.nFrames_sel)*100)) '% - ' num2str(toc) 's'])
        end

        temp.cdata          = read(vidObj, vidInfo.frameNums_sel(f));      % read frame
        vidInfo.tStamp(f)   = (vidInfo.frameNums_sel(f)-1)*(1/vidInfo.frameRate);
        Im                  = temp.cdata;
        if f == configVar.startFrame | countF == 1
            vidInfo.resolution = size(Im,1:2);
            textSize = round(vidInfo.resolution(1)/50);
            pupilCenterRef = vidInfo.resolution/2; % looks for pupil in center of image... adjust if necessary
            imDiv = []; 
        elseif f < configVar.startFrame+4 % adjust pupilCenterRef a couple more times after first time pupil center is detected
            pupilCenterRef = pupilCenter;
            imDiv = [];
        elseif pupilFitPerf > configVar.maxCircleFitError % bad pupil fit or not enough points detected
            if nFramesNoPupilDetected>configVar.waitForNFramesAfterPupilMiss 
                pupilCenterRef = vidInfo.resolution/2;
            else
                pupilCenterRef = pupilCenter;
            end

            imDiv = [];

        end
        ImGray              = rgb2gray(Im);
        ImOri               = ImGray;

    %% pupil detection

        % fast and accurate (straylight) with dynamically adaptating pupil center
        [B,R,pupilCenter,Res,imDiv] = detectPupil(ImGray,imDiv,pupilCenterRef,configVar.pupilRadiusRange,configVar.starburst.nThBins,configVar.startburst.stepSize);

        eyeData(f,1)             = sum(Im(:)==255);
        eyeData(f,2)             = R;
        eyeData(f,3)             = pupilCenter(2);
        eyeData(f,4)             = pupilCenter(1);
        eyeData(f,5)             = Res;
        eyeData(f,6)             = vidInfo.tStamp(f);
        if length(B) > 0
            eyeData(f,7)             = size(B,1);
        end


    %% Track SURF features and store their radius

        pupilFitPerf = eyeData(f,5);

        % redetection of track points happens when not enough
        % points are visible or when the pupil is not detected
        % properly. Only when the pupil is not detected properly (often due to a blink), a
        % couple frames is waited until tracking points are
        % detected again because eyelid needs to go back up first 

        % redetect features ... only when pupil is detected (e.g., after a blink)
        % and only when 10 frames are ok.

        if pupilDetected & nFramesSeqPupilDetected==configVar.waitForNFramesAfterPupilMiss | pupilDetected & nFramesSeqPupilDetected>configVar.waitForNFramesAfterPupilMiss & pointsDetected==0 % pupil is detected for X frames
            nTimesPointsRedetected = nTimesPointsRedetected+1;
            points = detectMinEigenFeatures(ImOri);
            points = points.Location;
            if size(points,1) > configVar.nFeaturePointRange(2)
                selPointsIdx = randperm(size(points,1));
                points = points(selPointsIdx(1:configVar.nFeaturePointRange(2)),:);
            end
            initialize(pointTracker, points, ImOri);
            oldPoints = points;
            pointsDetected = 1;
        end


        if pointsDetected & pupilDetected % check change in location of points
            [points, isFound]   = step(pointTracker, ImOri);
            visiblePoints       = points(isFound, :);
            oldInliers          = oldPoints(isFound, :);

            numPts = size(visiblePoints,1);
        end

        if numPts < configVar.nFeaturePointRange(1) % not enough points visible in image
            points = [];
            oldPoints = [];
            release(pointTracker);
            pointsDetected = 0;
        end

        if pupilFitPerf > configVar.maxCircleFitError % no good pupil fit
            points = [];
            oldPoints = [];
            release(pointTracker);
            pointsDetected = 0;

            pupilDetected = 0;
            nFramesSeqPupilDetected = 0;
            nFramesNoPupilDetected = nFramesNoPupilDetected+1;
        elseif pupilFitPerf <= configVar.maxCircleFitError  % correct pupil detection
            pupilDetected = 1;
            nFramesSeqPupilDetected = nFramesSeqPupilDetected+1;
            nFramesNoPupilDetected = 0;
        end


        if f == round(vidInfo.nFrames_sel) % end of video; reset, otherwise pointtracker cannot be initialized for next video
            release(pointTracker);
        end


    %% Calculate ratio of differences in radius of SURF points between inner and outer iris annulus 

        if numPts > configVar.nFeaturePointRange(1) & pupilFitPerf <= configVar.maxCircleFitError

            % loop through radii (inner vs outer circle)
            % and calculate difference in radius between all points
            % within inner and outer circle, and calc ratio of average
            % distance
            visiblePointsNorm = visiblePoints-repmat([pupilCenter(2) pupilCenter(1)],numPts,1); % with respect to pupil center

            [visiblePointsNormTH,visiblePointsNormR] = cart2pol(visiblePointsNorm(:,1),visiblePointsNorm(:,2));

            % only further than the pupil border
            tempR = visiblePointsNormR;
            selectInsideIrisVect = tempR>eyeData(f,2) & tempR<configVar.scleraRadius;

            visiblePointsNormR = visiblePointsNormR(selectInsideIrisVect);
            visiblePointsNormTH = visiblePointsNormTH(selectInsideIrisVect);
            visiblePointsNorm = visiblePointsNorm(selectInsideIrisVect,:);
            visiblePointsOri = visiblePoints(selectInsideIrisVect,:); % for stability calculations
            nPointsInIris = size(visiblePointsNorm,1);
            
            % average radial distance between all iris features
            D = abs(bsxfun(@minus,visiblePointsNormR, visiblePointsNormR.')); % square form
            aveD = mean(D,2);
            eyeData(f,9)             = mean(aveD);


            [ZTemp, RTemp, ~] = fitcircle(visiblePointsOri);
            eyeData(f,30) = ZTemp(1); %store average x position
            eyeData(f,43) = ZTemp(2); %store average y position

            % calculate feature distance within 10 iris annuli, each located at a different radius
            % also calculate center based on features per annulus
            radSelPrc = linspace(0,100,11); % radius divisions based on percentile ... not correct
            for i = 1:10
                % calculate radial (eccentricity) distance between points
                dData = visiblePointsNormR(visiblePointsNormR>=prctile(visiblePointsNormR,radSelPrc(i)) & visiblePointsNormR<prctile(visiblePointsNormR,radSelPrc(i+1)));
                D0 = abs(bsxfun(@minus,dData, dData.')); % square form
                aveD0 = mean(D0,2);
                eyeData(f,9+i) = mean(aveD0);
                eyeData(f,19+i) = length(dData); %store number of points tracked per annulus

                % fit circle and calculate x-y coordinates of circle
                [ZTemp, RTemp, ~] = fitcircle(visiblePointsOri(visiblePointsNormR>=prctile(visiblePointsNormR,radSelPrc(i)) & visiblePointsNormR<prctile(visiblePointsNormR,radSelPrc(i+1)),:));
                eyeData(f,32+i) = ZTemp(1); %store average x position
                eyeData(f,45+i) = ZTemp(2); %store average y position
                
                        
            end

            [ZTemp, RTemp, ~] = fitcircle(visiblePointsOri(visiblePointsNormR<prctile(visiblePointsNormR,50),:));
            eyeData(f,31) = ZTemp(1); %store average x position of pupillary, inner features
            eyeData(f,44) = ZTemp(2); %store average y position of pupillary, inner features
            
            [ZTemp, RTemp, ~] = fitcircle(visiblePointsOri(visiblePointsNormR>=prctile(visiblePointsNormR,50),:));
            eyeData(f,32) = ZTemp(1); %store average x position of ciliary, outer features
            eyeData(f,45) = ZTemp(2); %store average y position of ciliary, outer features

        else
            nPointsInIris = 0;
            eyeData(f,9:end) = NaN;
        end

        eyeData(f,8) = nPointsInIris;

        %% Show image processing result

        if configVar.disp.saveVideoIrissometryOutput | configVar.disp.inspectFramesManualPace | configVar.disp.inspectFrames
            if size(B,1) > 4
                % Insert the bounding box around the object being tracked
                bTemp = fliplr(B)';

                % highlight pupil border
                ImOri = insertShape(ImOri, 'Polygon', bTemp(:)','Opacity',0.8,'Color','red');
                ImOri = insertMarker(ImOri, bTemp', '+','Color', 'red');

                % highlight pupil area
                ImOri = insertShape(ImOri, 'FilledPolygon', bTemp(:)','Opacity',0.2,'Color','red');

                % show pupil circle fit
                ImOri = insertShape(ImOri, 'Circle', [pupilCenter(2) pupilCenter(1) R],'Opacity',0.8,'Color','green');

                % show tracked feature points 
                if numPts > configVar.nFeaturePointRange(1) & pupilFitPerf < configVar.maxCircleFitError  
                    ImOri = insertMarker(ImOri, visiblePointsNorm+repmat([pupilCenter(2) pupilCenter(1)],nPointsInIris,1), '+','Color', 'green','Size',5);
                end
            else
                ImOri = insertText(ImOri,[pupilCenter(2)*0.15 pupilCenter(1)*0.1],{'No pupil detected! Maybe adjust pupilRadiusRange in config variable'}, 'FontSize', textSize*2,'BoxColor','red','TextColor','w');
            end
            
            % show preset limbus radius and pupil search range
            ImOri = insertShape(ImOri, 'Circle', [pupilCenter(2) pupilCenter(1) configVar.scleraRadius],'Opacity',0.8,'Color','blue');
            ImOri = insertText(ImOri,[pupilCenter(2) pupilCenter(1)+configVar.scleraRadius],{'Limbus'}, 'FontSize', textSize,'BoxColor','blue','TextColor','w');
            
            ImOri = insertShape(ImOri, 'Circle', [pupilCenter(2) pupilCenter(1) configVar.pupilRadiusRange(1)],'Opacity',0.8,'Color','cyan');
            ImOri = insertText(ImOri,[pupilCenter(2) pupilCenter(1)+configVar.pupilRadiusRange(1)],{'Min pupil'}, 'FontSize', textSize,'BoxColor','cyan');
            
            ImOri = insertShape(ImOri, 'Circle', [pupilCenter(2) pupilCenter(1) configVar.pupilRadiusRange(2)],'Opacity',0.8,'Color','cyan');
            ImOri = insertText(ImOri,[pupilCenter(2) pupilCenter(1)+configVar.pupilRadiusRange(2)],{'Max pupil'}, 'FontSize', textSize,'BoxColor','cyan');


            % show center of pupil
            ImOri = insertMarker(ImOri, fliplr(pupilCenter), '+','Color', 'red');

            textString = {['Fit res: ' num2str(pupilFitPerf)],...
                          ['Frame #: ' num2str(f)],...
                          ['# Frames pupil: ' num2str(nFramesSeqPupilDetected)],...
                          ['# points image: ' num2str(numPts)],...
                          ['# points iris: ' num2str(nPointsInIris)],...
                          ['Time: ' num2str(round(vidInfo.tStamp(f),2))]};
                      
            textPosition = zeros(size(textString,2),2);
            for i = 1:size(textString,2)
%                 textPosition(i,:) = [vidInfo.resolution(1)*textSize vidInfo.resolution(2)*textSize+2*vidInfo.resolution(2)*textSize*(i-1)];
                textPosition(i,:) = [10 (textSize*2)*i];
            end
            
            ImOri = insertText(ImOri,textPosition,textString, 'FontSize', textSize);

        end
        if configVar.disp.saveVideoIrissometryOutput
            writeVideo(vidObj_fit,uint8(ImOri))
        end

        if  ( configVar.disp.inspectFrames )
            figure(1001);
            if f == configVar.startFrame
                set(gcf,'Position',[0 0 1920 1080])
            end
            imshow(ImOri);

            if configVar.disp.inspectFramesManualPace
                WaitForKeyPress({'SPACE'});
            end
        end
    end
    if configVar.disp.saveVideoIrissometryOutput
        close(vidObj_fit)
    end
    save([vidFileName '_output.mat'],'eyeData','eyeDataHeader','vidInfo');
    disp(['Saved output to: ' vidFileName '_output.mat'])
else
    disp('Analysis output file already exists.')
    disp(['Loading: ' vidFileName '_output.mat'])
    load([vidFileName '_output.mat'],'eyeData','eyeDataHeader','vidInfo');
end
