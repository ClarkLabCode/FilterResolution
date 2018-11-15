%% script shows intuition for super-resolution
clear
clc
rng(2); % to get similar behavior each time

sampleInterval = 50; % sample every 500 ms

t = [0:500];
B = (1-exp(-t/2)).*(exp(-t/10)/10 - t.*exp(-t/20)/20/20);
A = 1;

% How far forward and back to calculate filter, in ms
extentForwards = 150;
extentBackwards = 10;

%% Show how more samples is better
impulse = zeros(extentBackwards+1+extentForwards,1);
impulse(extentBackwards+1) = 1;
impulseResponse = filter(B,A,impulse);

%Run through the same number of samples to generate exactly the same noise
%as figure3
for tLength = [1e5 5e4 1.5e4]
    tVec = 1:tLength;
    
    stim = randn(size(tVec));
    resp = filter(B,A,stim);
    resp = resp + std(resp)*randn(size(tVec)); % this sets SNR = 1!
end

sampleIdxs = 1:sampleInterval:max(tVec);

isSampled = zeros(size(resp));
isSampled(sampleIdxs) = 1;

responseSampled = resp.*isSampled;

responseSampled = responseSampled - mean(responseSampled);
stim = stim - mean(stim);


isSampled = zeros(size(resp));
isSampled(sampleIdxs) = 1;
stimMatrixForInterp = toeplitz(stim',zeros(1,extentForwards+extentBackwards+1));
stimMatrixForInterp = stimMatrixForInterp(extentForwards+extentBackwards+1:end,:);

stimMatrixForSuper = toeplitz(stim',zeros(1,extentForwards+extentBackwards+1+sampleInterval*2));
stimMatrixForSuper = stimMatrixForSuper(extentForwards+extentBackwards+1+sampleInterval*2:end,:);

figure(1);
clf;
responseInterped = interp1(sampleIdxs,responseSampled(logical(isSampled)),tVec);

numSamples = sum(isSampled(extentForwards+1+sampleInterval:end-extentBackwards-sampleInterval));
normalXcorr = responseSampled(1,extentForwards+1+sampleInterval:end-extentBackwards-sampleInterval)*stimMatrixForSuper/numSamples;
sampFreq = 1/sampleInterval;
triangleFilter = [0:sampFreq:1, 1-sampFreq:-sampFreq:0];
triangleFilter = triangleFilter./sum(triangleFilter);
normalXcorrFiltered = conv(normalXcorr,triangleFilter,'valid');

numSamplesInterp = length(isSampled(extentForwards+1+sampleInterval:end-extentBackwards-sampleInterval));
interpedXcorr = responseInterped(1,extentForwards+1+sampleInterval:end-extentBackwards-sampleInterval)*stimMatrixForSuper/numSamplesInterp;
deltaFilter = zeros(size(triangleFilter));
deltaFilter(ceil(length(deltaFilter)/2)) = 1;
interpedXcorr = conv(interpedXcorr,deltaFilter,'valid');

hold on;
plot([-extentBackwards:extentForwards]*10,[impulseResponse,normalXcorr(sampleInterval+1:end-sampleInterval)'])
plot([-extentBackwards:extentForwards]*10,interpedXcorr','LineWidth',2);
plot([-extentBackwards:extentForwards]*10,normalXcorrFiltered');
hold off
title('Triangle Filtering VS Interpolation XCorr')
legend({'Ground Truth','No Filtering','Interpolated','Super Res Triangle Filtered'})
set(gca,'xlim',[-extentBackwards,extentForwards]*10);