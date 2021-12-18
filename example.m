clear all
close all
clc


%% README
%
% This file contains an example of how to run the main irissometry function
% to analyze a close-up video of an eye and track feature within the iris
% frame to frame.
%
% This example file opens a pop-up window to select a directory with video 
% files to be analyzed one by one.
% 
% To understand what the input and output of the function is, enter "help
% irissometry" in the command.
% 
% REQUIRED MATLAB TOOLBOXES: 
% - computer vision toolbox
% - image processing toolbox
% - signal processing toolbox
% - statistics and machine learning toolbox
%
% Contact: marnixnaber@gmail.com

%% Planned improvements to the code

% - Automated limbus detection
% - Improved pupil border detection (using sophisticated computer vision techniques)

%% Set configuration parameters

% enter "help irissometry" in command for configuration details

configVar = struct();
configVar.scleraRadius = 275;                   % location of limbus, set manually [pixels]
configVar.pupilRadiusRange = [50 150];          % search range for pupil border [minimum maximum] in [pixels]
configVar.startFrame = 1;                       % first video frame to be analyzed
configVar.maxCircleFitError = 0.16;             % maximum allowed error in circle fitted to pupil border (otherwise blink is detected)
configVar.nFeaturePointRange = [600 2000];      % range of number of feature to be detected in iris [minimum maximum] 
configVar.waitForNFramesAfterPupilMiss = 10;    % number of frames
configVar.starburst.nThBins = 16;               % number of points detected on pupil border
configVar.startburst.stepSize = 5;              % distance in [pixels] to compare radial luminance increases to detect pupil border
configVar.disp.inspectFrames = 1;               % set to 1 to display popup with irissometry output
configVar.disp.inspectFramesManualPace = 0;     % set to 1 to display popup with irissometry output per frame, use space to show a subsequent frame
configVar.disp.saveVideoIrissometryOutput = 1;  % set to 1 to save a new video with irissometry output
configVar.overwriteOutputData = 0;              % set to 1 to overwrite existing irissometry output (_output.mat files)

%% select folder with video files

selectFileBy = 'directory'; % Options: 'directory' or 'files'
selectFileBy = 'files'; 
videoExtension = '*.mp4';
    
if strcmp(selectFileBy,'directory')
    videosFolder = uigetdir;
    tempFileStruct = dir(fullfile(videosFolder, [videoExtension]));
    videoFileNames = {};
    for k = 1:length(tempFileStruct)
        pathName = [tempFileStruct.folder '\'];
        videoFileNames{k} = tempFileStruct.name;
    end
    
elseif strcmp(selectFileBy,'files')
    [videoFileNames, pathName] = uigetfile(videoExtension,'Select the INPUT DATA FILE(s)','MultiSelect','on');
    if ~iscell(videoFileNames)
        videoFileNames = {videoFileNames};
    end
else
    disp('ERROR: string in variable selectFileBy not recognized. Should be "directory" or "files".')
end

%% run irissometry per video

for vidNum = 1:length(videoFileNames)
    % Automated limbus detection
    % WILL BE IMPLEMENTED LATER
    
    vidFileName = [pathName videoFileNames{vidNum}];
    [eyeData,eyeDataHeader,vidInfo] = irissometry(vidFileName,configVar);
    
    % output data in csv file
    T = array2table(eyeData, 'VariableNames',[eyeDataHeader]);
    T.Properties.VariableNames = eyeDataHeader;
    writetable(T,[vidFileName '_output.csv'],'Delimiter',',') 
    disp(['Saved output to: ' vidFileName '_output.csv'])
    
    % remove no iris detection episodes
    selectBool = eyeData(:,8)==0;
    eyeData(selectBool,:) = NaN;
    
    
    % Blink epoch removal and interpolation
    % WILL BE IMPLEMENTED LATER
    
    
    % plot some results
    plotResults;
    
end