% script to extract kernels from ArcLD recordings
clear all;

filterLength = 36; % in frames; 300 ms in our data
numFramesForward = 6; % in frames; 50 ms in our data

% Plots up 4 analyses, depending on choosePlot below
choosePlot = 4; % 1-4, depending, see below

switch choosePlot
    case 1 % Line scans (figure 4 supplement 1)
        % Load Data
        D = load('ArcLightLineScanData.mat');
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        pixelsPerMicron = D.pixelsPerMicron;
        oneSecondSegment = D.oneSecondSegment;
        
    case 2 % 2d data, treating each line as a sample (figure 4c)
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
           % and interpolating between frame-to-frame samples
           % (figure 4d pink)
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
        
    case 4 % 2d data, as 3, but supposing frames taken at 1/n original
           % frame rate (but equal frame duration) (figure 4d green)
        % Load Data
        
        D = load('ArcLight2dScanROI.mat');
        
        allResp = D.allRespForKernel;
        allStim = D.allStimForKernel;
        meanImage = D.meanMovie;
        roiMask = D.roiMask;
        pixelsPerMicron = D.pixelsPerMicron;
        
        sampleN = 13; % choose every N samples of the response to 
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

if choosePlot ~= 2
    for ii=1:size(allResp,2)
        % We'll be using bootstrapping to calculate accurate error bars for
        % the first ROI. (These are the error bars plotted in figure 4d)
        bootStrap = (ii == 1) & choosePlot ~= 1;
        [kernels(:,ii), errors(:,ii)] = extractKernel(allStim, allResp(:,ii), filterLength, numFramesForward,bootStrap);
    end
else
    % In the by-line case, we have to modify the extraction process to
    % correctly combine each line of the ROI
    for rr=1:length(roiIdx)
        f = find(roiMember == roiIdx(rr));
        % We'll be concatenating the stimulus and response matrix for each
        % line belonging to a particular ROI
        RBig = [];
        SBig = [];
        for ff = f
            stimulus = allStim;
            stimulus = stimulus-nanmean(stimulus); % mean subtract it!
            response = allResp(:,ff);

            % Choose which response indices to analyze, excluding all nans
            chooseResponseInds = find(~isnan(response));
            chooseResponseInds = chooseResponseInds(chooseResponseInds>filterLength);
            chooseResponseInds = chooseResponseInds(chooseResponseInds<length(stimulus)-numFramesForward);

            % Construct stimulus and response matrices from selected response
            % indices
            S = zeros(length(chooseResponseInds), filterLength+numFramesForward+1);
            R = zeros(length(chooseResponseInds), 1);

            for ii=1:length(chooseResponseInds)
                ch = chooseResponseInds(ii); % current frame
                R(ii) = response(ch);
                S(ii,:) = stimulus(ch-filterLength:ch+numFramesForward);
            end
            RBig = [RBig;R];
            SBig = [SBig;S];
        end

        RBig = RBig-mean(RBig);

        % Use bootstrap method to find error bars for the first ROI
        bootStrap = (ii == 1);
        if bootStrap
            kernel = regress(RBig,SBig);
            CIs = bootci(100,{@regress,RBig,SBig},'alpha',1-0.682,'type',cper);
        else
            [kernel, CIs] = regress(RBig,SBig,1-0.682);
        end
        error = CIs(:,2)-kernel;

        % Reverse kernels for standard notation
        kernels(:,rr) = kernel(end:-1:1);
        errors(:,rr) = error(end:-1:1);
    end
end

%% Plot data for all cases

figure; hold on;
subplot(2, 1, 1);
hold on;
ts = [-numFramesForward:filterLength]'*1000/120;
ks = -kernels*(120);
es = -errors*(120);
PlotXvsY(ts,ks,'error',es);
plot([-numFramesForward,filterLength]/120,[0 0],'k:');
set(gca,'ylim',[-4 10]);
xlabel('time (ms)');
ylabel('-dF/F/contrast/s');
title(['All kernels; method ' num2str(choosePlot)]);

% draws a line indicating the sampling interval
plot([0 1000/13],-1*[1 1],'k-','linewidth',2);

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