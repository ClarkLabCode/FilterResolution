%% script shows intuition for super-resolution
clear;
rng(2); % to get similar behavior each time

sampleInterval = 50; % sample every 500 ms

t = [0:500];
B = (1-exp(-t/2)).*(exp(-t/10)/10 - t.*exp(-t/20)/20/20);
A = 1;

% How far forward and back to calculate filter, in ms
extentForwards = 150;
extentBackwards = 10;

%% Show how more samples is better
figure(1);
impulse = zeros(extentBackwards+1+extentForwards,1);
impulse(extentBackwards+1) = 1;
impulseResponse = filter(B,A,impulse);

lengthsVec = 10.^(4:(1/6):7);
for ii = 1:length(lengthsVec)
    tVec = 1:lengthsVec(ii);
    
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
    
    error(ii) = sqrt(mean((sampledKernel-impulseResponse).^2))./sqrt(mean(impulseResponse.^2));
    numSamples(ii) = length(validSampledResp);
end
title('Noise Scaling with # of samples');
loglog(numSamples,error,...
       numSamples,15./sqrt(numSamples));
xlabel('number of samples')
ylabel('relative noise magnitude')
colors = get(gca, 'ColorOrder');
text(4.0e4,1/10,'${1}\over{\sqrt{N}}$','Color',colors(2,:),'Interpreter','latex','FontSize',14)

function toeMatrix = calcSparseToeplitz(stim,sampleIdxs,extentForwards,extentBackwards)
sampleTimesSel = sampleIdxs( (sampleIdxs > extentForwards) & (sampleIdxs < length(stim)-extentBackwards+1));
toeMatrix = zeros(length(sampleTimesSel),extentBackwards+extentForwards+1);
for ii=1:length(sampleTimesSel)
    % note, this is just an explicit computation of the cross-correlation
    toeMatrix(ii,:) = stim(sampleTimesSel(ii)+[-extentForwards:extentBackwards]);
end
toeMatrix = fliplr(toeMatrix);
end
