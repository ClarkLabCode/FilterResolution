% Create a demo impulse response function to use with MATLAB's "filter"
% command.
t = [0:500];
B = (1-exp(-t/2)).*(exp(-t/10)/10 - t.*exp(-t/20)/20/20);
A = 1;

% set how far forward and back to calculate filter, in number of samples
extentPast = 150;
extentFuture = 10;

stimEnd = 1e4; % Length of experiment in seconds
stimStart = 15; % The time when the stimulus turns on
stimFrequency = 120; % Stimulus update rate

numRespSamples = 10000;

% Define when the stimulus is presented and response is sampled
stim_ts = stimStart:1/stimFrequency:stimEnd;
% Random response sampling:
% resp_ts = stimEnd*rand(numRespSamples,1);
% Uniform response sampling:
resp_ts = linspace(0,stimEnd,numRespSamples);

% Calculate ground truth response
impulse = zeros(extentFuture+1+extentPast,1);
impulse(extentFuture+1) = 1;
impulseResponse = filter(B,A,impulse);

% Generate high resolution stimulus and response
stim = randn(size(stim_ts));
stimFull = [zeros(1,stimStart*stimFrequency), stim];
respFull_ts = 0:1/stimFrequency:stimEnd;
respFull = filter(B,A,stimFull);
respFull = respFull + std(respFull)*randn(size(respFull)); % this sets SNR = 1!
resp = interp1(respFull_ts,respFull,resp_ts);

%% Demo!
[f_hat,f_hat_ts] = calculateVoxelTimingFilter(stim_ts,stim,resp_ts,resp,extentPast,extentFuture,'ols');

plot(f_hat_ts,f_hat,f_hat_ts,impulseResponse)
legend({'$\hat{f}$','ground truth'},'Interpreter','latex')

