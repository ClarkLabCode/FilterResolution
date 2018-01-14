function [kernel, error] = extractKernel(stim,resp,numFramesBack,numFramesForward,bootStrap)

% Function takes arguments
% stim: column vector of stimulus, at a high frame rate, typically
% resp: response for each sample of stim; nan whanere
% there was no sampling
% numFramesBack: number of frames back in time to compute the filter
% numFramesForward: number of frames forward in time to compute the filter
% Output:
% kernel: the number of elements is (numFramesBack+numFramesForward+1), since
% the response with no delay is not counted among Forward or Back frames
% here

% Assign stimulus
stimulus = stim;
stimulus = stimulus-nanmean(stimulus); % mean subtract it!
response = resp;

% Choose which response indices to analyze, excluding all nans
chooseResponseInds = find(~isnan(response));
chooseResponseInds = chooseResponseInds(chooseResponseInds>numFramesBack);
chooseResponseInds = chooseResponseInds(chooseResponseInds<length(stimulus)-numFramesForward);

% Construct stimulus and response matrices from selected response
% indices
S = zeros(length(chooseResponseInds), numFramesBack+numFramesForward+1);
R = zeros(length(chooseResponseInds), 1);

for ii=1:length(chooseResponseInds)
    ch = chooseResponseInds(ii); % current frame
    R(ii) = response(ch);
    S(ii,:) = stimulus(ch-numFramesBack:ch+numFramesForward);
end

R = R-mean(R);
if bootStrap
    kernel = regress(R,S);
    CIs = bootci(100,{@regress,R,S},'alpha',1-0.682,'type','cper')';
else
    [kernel, CIs] = regress(R,S,1-0.682);
end
error = CIs(:,2)-kernel;
% Could compute Rhat here easily by finding Rhat = S*k;

kernel = kernel(end:-1:1); % Reverse them for standard notation
error = error(end:-1:1);
end
