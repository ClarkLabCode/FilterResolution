%% script shows intuition for super-resolution

rng(2); % to get similar behavior each time

sampleInterval = 50; % sample every 500 ms

t = [0:500];
B = (1-exp(-t/2)).*(exp(-t/10)/10 - t.*exp(-t/20)/20/20);
A = 1;

% How far forward and back to calculate filter, in ms
extentForwards = 150;
extentBackwards = 10;

%% Show how more samples is better
figure(2);
impulse = zeros(extentBackwards+1+extentForwards,1);
impulse(extentBackwards+1) = 1;
impulseResponse = filter(B,A,impulse);

plot([-extentBackwards:extentForwards]*10,impulseResponse,'color',0.5*[1 1 1],'linew',2);
hold all;
for tLength = [1e5 5e4 1.5e4]
    tVec = 1:tLength;
    
    stim = randn(size(tVec));
    resp = filter(B,A,stim);
    resp = resp + std(resp)*randn(size(tVec)); % this sets SNR = 1!
    
    sampleIdxs = 1:sampleInterval:max(tVec);
    
    isSampled = zeros(size(resp));
    isSampled(sampleIdxs) = 1;
    
    responseSampled = resp.*isSampled;
    
    validSampleIdxs = sampleIdxs( (sampleIdxs > extentForwards) & (sampleIdxs < length(resp)-extentBackwards+1));
    validSampledResp = resp(validSampleIdxs);
    validSampledResp = validSampledResp - mean(validSampledResp);
    stimMatrix = calcSparseToeplitz(stim,sampleIdxs,extentForwards,extentBackwards);
    stimMatrix = stimMatrix-mean(stimMatrix);
    sampledKernel = regress(validSampledResp',stimMatrix);
    plot([-extentBackwards:extentForwards]*10,sampledKernel);
end
hold off;
title('Noise Scaling with # of samples');
legend({'Ground Truth','1e5 samples','5e4 samples','1.5e4 samples'})
set(gca,'xlim',[-extentBackwards,extentForwards]*10);

%% Show simulated stimulus and response
figure (1)
plotRange = 125:325;
subplot(3,1,1)
plot(tVec(plotRange)*10,stim(plotRange))
axis('tight')
ylabel('stimulus');
subplot(3,1,2);
cla
hold on;
for sampleIdx = find(isSampled(plotRange))
    tVecPlot = tVec(plotRange)*10;
    plot(tVecPlot(sampleIdx)*[1 1],[0 1],'k');
end
hold off;
xlim(10*[min(tVec(plotRange)),max(tVec(plotRange))]);
ylim([0 1]);
ylabel('response samples');
subplot(3,1,3);
isSampledInRange = zeros(size(tVec));
isSampledInRange(plotRange) = isSampled(plotRange);
isSampledInRange = logical(isSampledInRange);
plot(tVec(plotRange)*10,resp(plotRange),tVec(isSampledInRange)*10,resp(isSampledInRange),'o');
axis('tight')
ylabel('response');
xlabel('time (ms)');
%% Show how gaussian filtering improves your estimate
figure(3)
plot([-extentBackwards:extentForwards]*10,impulseResponse,'color',0.5*[1 1 1],'linew',2);
hold all;
for filterWidth = [1 10 25] %[1 5 10 15 20]
    plot([-extentBackwards:extentForwards]*10,smoothdata(sampledKernel,'gaussian',filterWidth));
    % note: Gaussian has a std of window/5, found empirically because
    % Matlab documentation handily does not say
end
hold off;
title('Gaussian Filtering Reduces Noise');
legend({'Ground Truth','No Filtering','20 ms std Gaussian filter','50 ms std Gaussian filter'});
set(gca,'xlim',[-extentBackwards,extentForwards]*10);

%% Show how triangle filtering is equivalent to interpolating
figure(4);

isSampled = zeros(size(resp));
isSampled(sampleIdxs) = 1;
stimMatrixFull = toeplitz(stim',zeros(1,extentForwards+extentBackwards+1));
stimMatrixFull = stimMatrixFull(extentForwards+extentBackwards+1:end,:);


responseInterped = interp1(sampleIdxs,resp(logical(isSampled)),tVec);
numSamplesInterp = length(isSampled(extentForwards+1:end-extentBackwards));
interpedKernel = regress(responseInterped(1,extentForwards+1:end-extentBackwards)',stimMatrixFull);
sampFreq = 1/sampleInterval;
triangleFilter = [0:sampFreq:1, 1-sampFreq:-sampFreq:0];
triangleFilter = triangleFilter./sum(triangleFilter);
sampledKernelFiltered = conv(sampledKernel,triangleFilter,'same');

plot([-extentBackwards:extentForwards]*10,[impulseResponse,sampledKernel,interpedKernel,sampledKernelFiltered]);
title('Triangle Filtering VS Interpolation OLS')
legend({'Ground Truth','No Filtering','Interpolated','Super Res Triangle Filtered'})
set(gca,'xlim',[-extentBackwards,extentForwards]*10);


%%
figure(5);

% Laguerre
laguerreFuncs = getLaguerrePolys(length(sampledKernel(extentBackwards+1:end)),5,0.8);
laguerreCoeffs = sampledKernel(extentBackwards+1:end)'*laguerreFuncs;
outputKernel = laguerreFuncs*laguerreCoeffs';

laguerreEst = zeros(size(sampledKernel));
laguerreEst(extentBackwards+1:end) = outputKernel;

% ASD
thisPath = fileparts(mfilename('fullpath'));
addpath(genpath(thisPath));
minlen=0;

[asdEstFull,~] = fastASD(stimMatrix,validSampledResp',length(sampledKernel),minlen);

% Plot
plot([-extentBackwards:extentForwards]*10,[impulseResponse,sampledKernel,laguerreEst,asdEstFull]);
title('Regularization Preserves Enhances SNR w/o sacrificing Resolution')
legend({'Ground Truth','No Regularization','Laguerre Poly Regularization','ASD Regularization'});
set(gca,'xlim',[-extentBackwards,extentForwards]*10);

function toeMatrix = calcSparseToeplitz(stim,sampleIdxs,extentForwards,extentBackwards)
    sampleTimesSel = sampleIdxs( (sampleIdxs > extentForwards) & (sampleIdxs < length(stim)-extentBackwards+1));
    toeMatrix = zeros(length(sampleTimesSel),extentBackwards+extentForwards+1);
    for ii=1:length(sampleTimesSel)
        % note, this is just an explicit computation of the cross-correlation
        toeMatrix(ii,:) = stim(sampleTimesSel(ii)+[-extentForwards:extentBackwards]);
    end
    toeMatrix = fliplr(toeMatrix);
end
