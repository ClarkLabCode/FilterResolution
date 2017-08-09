% script to extract kernels from ArcLD recordings
clear all;

filterLength = 36; % in frames; 300 ms in our data
numFramesForward = 6; % in frames; 50 ms in our data

% Plots up 5 analyses, depending on choosePlot below
choosePlot = 2; % 1-5, depending, see below
switch choosePlot
    case 1 % Line scans
        % Load Data
        D = load('ArcLightLineScanData.mat');
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        pixelsPerMicron = D.pixelsPerMicron;
        oneSecondSegment = D.oneSecondSegment;
        
    case 2 % 2d data, treating each line as a sample
        % Load Data
        D = load('ArcLight2dScanByLine.mat');
        
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        origRoiMask = D.origRoiMask;
        pixelsPerMicron = D.pixelsPerMicron;
        
        % Associate each line scan segment with an original ROI
        roiIdx = unique(origRoiMask(origRoiMask>0));
        lineIdx = unique(roiMask(roiMask>0));
        for ii=1:length(lineIdx)
            f = find(roiMask == lineIdx(ii));
            roiMember(ii) = origRoiMask(f(1)); % assign it the index
        end
        
    case 3 % 2d data, treating each ROI as single sample at ROI mid point time
        % Load Data
        D = load('ArcLight2dScanROI.mat');
        
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        pixelsPerMicron = D.pixelsPerMicron;

        % allRespForKernel = stimLength x ROIs matrix -- response at ROI
        % mid point
        % allStimForKernel = stimLength x 1 vector
        % meanMovie = movieHeight x movieWidth matrix
        % roiMask = movieHeight x movieWidth matrix
        
    case 3.1 % 2d data, treating each ROI as single fast sample,
        % timing mapped to first point of frame
        
        D = load('ArcLight2dScanROIAlignToFirstFrame.mat');
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        pixelsPerMicron = D.pixelsPerMicron;
                
    case 4 % 2d data, as 3, but interpolating between frame-to-frame samples
        % Load Data
        D = load('ArcLight2dScanROI.mat');
        
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        pixelsPerMicron = D.pixelsPerMicron;
        
        for ii=1:size(allResp,2)
            allResp(:,ii) = interpolateNan(allResp(:,ii));
        end
        
    case 5 % 2d data, as 3, but supposing frames taken at 1/n original frame rate (but equal frame duration)
        % Load Data
        
        D = load('ArcLight2dScanROI.mat');
        
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        pixelsPerMicron = D.pixelsPerMicron;
        
        sampleN = 3; % choose every N samples of the response to 
        % simulate lower aquisition rate (but same linescan rate); for choosePlot 3
        for ii=1:size(allResp,2)
            locResp = allResp(:,ii);
            newResp = nan(size(locResp));
            f = find(~isnan(locResp));
            newResp(f(1:sampleN:end)) = locResp(f(1:sampleN:end));
            allResp(:,ii) = newResp; % replace with downsampled version
        end
        
end

%% Extract the kernels for all cases

for ii=1:size(allResp,2)
    kernels(:,ii) = extractKernel(allStim, allResp(:,ii), filterLength, numFramesForward);
end

if (choosePlot == 2)
    for ii=1:length(roiIdx) % average them by big ROI
        f = find(roiMember == roiIdx(ii));
        tempK(:,ii) = mean(kernels(:,f),2);
    end
    kernelsLine = kernels;
    kernels = tempK;
end

%% Plot data for all cases

figure; hold on;
subplot(2, 1, 1);
hold on;
plot([-numFramesForward:filterLength]/120*1000,-kernels*(120/1000));
plot([-numFramesForward,filterLength]/120*1000,[0 0],'k:');
set(gca,'ylim',[-2 6]*1e-3);
xlabel('time (ms)');
ylabel('-dF/F/contrast/ms');
title(['All kernels; method ' num2str(choosePlot)]);
switch choosePlot
    case {2 3 3.1 4}
        % draws a line indicating the sampling interval
        plot([0 1000/13],-0.001*[1 1],'k-','linewidth',2);
    case 5
        % draws a line indicating the sampling interval
        plot([0 1000/13*sampleN],-0.001*[1 1],'k-','linewidth',2);
end

% plot ROIs
subplot(2, 1, 2);
if choosePlot == 1
    meanImage = oneSecondSegment;
end

imagesc(meanImage);
axis off;axis tight; axis equal;colormap gray; hold on;

if size(roiMask, 1) == 1 %linescan
    roiMask = repmat(roiMask, size(oneSecondSegment, 1), 1);
end

roiMaskChoices = false(size(roiMask));
roiChoiceInds = unique(roiMask(:));
roiChoiceInds(roiChoiceInds==0) = [];
roiMaskOutlines = cell(0, 0);
roiInd = 1;

for i = 1:length(roiChoiceInds)
    roiMaskChoices(roiMask==roiChoiceInds(i)) = true;
    bndrs = bwboundaries(roiMaskChoices);
    roiMaskOutlines(end+1, 1:length(bndrs)) = bndrs;
    [indRows, indCols] = find(roiMaskChoices);
    roiCenterOfMass(i, 1) = mean(indRows);
    roiCenterOfMass(i, 2) = mean(indCols);
    roiMaskChoices(roiMask==roiChoiceInds(i)) = false;
    roiInd = roiInd+length(bndrs);
end
for i = 1:size(roiMaskOutlines, 1)
    for j = 1:size(roiMaskOutlines, 2)
        if ~isempty(roiMaskOutlines{i, j})
            lnOut = plot(roiMaskOutlines{i, j}(:, 2), roiMaskOutlines{i, j}(:, 1), 'LineWidth', 3);
        end
    end
    text(roiCenterOfMass(i, 2), roiCenterOfMass(i, 1), num2str(i), 'HorizontalAlignment', 'center', 'Color',  [1 1 1]);
end
% 10 micron scale bar
plot([10 10+10*pixelsPerMicron], [size(meanImage, 1)-10 size(meanImage, 1)-10], 'w', 'LineWidth', 4);