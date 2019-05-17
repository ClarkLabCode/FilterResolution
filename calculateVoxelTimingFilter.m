function [f_hat, f_hat_ts, S] = calculateVoxelTimingFilter(stim_ts,stim,resp_ts,resp,numStimPast,numStimFuture,method,minLen)
% CALCULATEVOXELTIMINGFILTER calculates an estimate of an impulse
% response function (a.k.a. filter or kernel) from a set of measurements.
% This estimate is super-resolved because the resolution of the estimate 
% depends only on the update rate of the stimulus, and not on the
% measurement frequency.
%
%   [f_hat, f_hat_ts] =
%   CALCULATEVOXELTIMINGFILTER(stim_ts,stim,resp_ts,resp,numStimPast,numStimFuture)
%   calculates this estimate using the default "cross-correlation" or
%   "reverse correlation" method.
%
%   INPUTS:
%
%   stim_ts and resp_ts are vectors that specify the time that each
%   stimulus and response was acquired. Responses are assumed to be
%   acquired in a very short amount of time. This function requires that
%   the stimuli are presented at regular intervals, but the responses
%   can be acquired at arbitrary times.
%
%   stim and resp are vectors that specify the values of the stimulus and
%   response at each element of stim_ts and resp_ts. stim and resp do not
%   need to be the same length as each other, but they do need to be the
%   same length as stim_ts and resp_ts respectively.
%
%   numStimPast is the number of elements of the stimulus presented before
%   the measurement of a response that are included in the filter estimate.
%   For instance, if numStimPast is equal to 3, this function will fit an
%   estimate of response kernel which can predict the effect of the three
%   preceeding stimuli plus the current stimulus on the currently measured
%   response.
%
%   Similarly, numStimFuture is the number of elements of the stimulus
%   presented after the measurement of the response that are used to
%   estimate that response. Usually a response cannot predict future
%   stimuli, so we expect that these elements of the kernel/filter/impulse
%   response function will be approximately zero.
%
%   OUTPUTS:
%
%   f_hat is the estimated impulse response function, deterimened at
%   f_hat_ts timepoints relative to the stimulus.
%
%   After calling this function, you can run plot(f_hat_ts,f_hat) to see
%   the resulting filter.
%
%
%   [f_hat, f_hat_ts] =
%   CALCULATEVOXELTIMINGFILTER(stim_ts,stim,resp_ts,resp,numStimPast,numStimFuture,method,minLen)
%   Also specifies a method, which can be 'xcorr' (cross-correlation or
%   reverse correlation) (DEFAULT), 'ols' (ordinary least squares), or 'asd'
%   (automatic smoothness determination). 'ols' takes longer to compute
%   than 'xcorr' and requires a well-conditioned stimulus (all frequencies
%   are represented), but can yield a less noisy estimate with fewer
%   samples and account for stimulus autocorrelation (stimuli that are not
%   "white"). 'asd' uses Automatic Smoothness Determination (Aoi & Pillow
%   2017). This requires installation of the fastASD function from their
%   <a href="https://github.com/pillowlab/fastASD">github page</a>. Make
%   sure it is in your MATLAB path.
%
%   If 'asd' is selected, the minLen parameter for the fastASD method can
%   be specified. The default is 0.
%
%
%   [f_hat, f_hat_ts, S] =
%   CALCULATEVOXELTIMINGFILTER(stim_ts,stim,resp_ts,resp,numStimPast,numStimFuture,...)
%   Also outputs a stimulus matrix, S, such that f_hat'*S is a
%   reconstruction of the response.

    %% Input sanitation

    if nargin < 7
        method = 'ols';
    end
    
    if nargin < 8
        minLen = 0;
    end
    
    assert(isvector(stim_ts) && isvector(stim) && isvector(resp_ts) && isvector(resp),...
           'stim_ts, stim, resp_ts, and resp must all be vectors');
    % Make all inputs column vectors
    stim_ts = stim_ts(:);
    stim = stim(:);
    resp_ts = resp_ts(:);
    resp = resp(:);
    
    assert(length(stim_ts) == length(stim),['stim_ts and stim inputs must' ...
                                            'be of equal length']);
    assert(length(resp_ts) == length(resp),['resp_ts and resp inputs must' ...
                                            'be of equal length']);
    
    assert(all(diff(stim_ts) > 0),'stim_ts must be monotonically increasing');
                                        
    assert(numStimPast >=0,'numStimBack must be nonnegative');
    assert(numStimFuture >=0,'numStimForward must be nonnegative');
    
    if ~ischar(method) || ~any(strcmpi(method,{'xcorr','ols','asd'}))
        error('please enter a valid method string: ''xcorr'',''ols'', or ''asd''');
    end
    method = lower(method);
    
    if strcmp(method,'asd') && ~exist('fastASD','file')
        error('could not find fastASD function. Please make sure it has been added to your path');
    end
    
    if ~strcmp(method,'asd') && nargin > 7
        warning('minlen argument is not used when method is not ''asd''');
    end
    
    stimDt = mean(diff(stim_ts));
    stimDtStd = std(diff(stim_ts));
    assert(stimDtStd/stimDt < 0.1,['stimuli must be presented at regular'...
                                   'intervals for this analysis.' ...
                                   'Deviations from this will cause' ...
                                   'distortions in the extrated filter']);
    
    %% Create S and r
    
    % get a more accurate interpolation of stimulus timing
    b = [(1:length(stim_ts))',ones(length(stim_ts),1)]\stim_ts;
    timeOffset = b(2);
    stimDt = b(1);

    % find the stim idxs that correspond to the response times. Assume that
    % the stimulus is presented at a contstant rate
    
    sIdxs =  round((resp_ts-timeOffset)/stimDt);
    validIdxs = (sIdxs-numStimPast > 0) & (sIdxs+numStimFuture <= length(stim));
    sIdxs = sIdxs(validIdxs);
    
    fastestMethod = find(min(predictSTime((numStimPast+numStimFuture+1),length(resp),length(stim))));
    switch fastestMethod
        case 1
            S = generateS_sparse(stim,sIdxs,numStimFuture,numStimPast);
        case 2
            S = generateS_sparseTranspose(stim,sIdxs,numStimFuture,numStimPast);
        case 3
            S = generateS_dense(stim,sIdxs,numStimFuture,numStimPast);
    end
    
    r = resp(validIdxs);
    
    %% Calculate f_hat
    
    switch method
        case 'xcorr'
            f_hat = r' * S / length(r);
        case 'ols'
            f_hat = S\r;
        case 'asd'
            f_hat = fastASD(S,r,(numStimFuture+numStimPast+1),minLen);
    end
    
    if fastestMethod < 3 % See comment above generateS_dense
        f_hat = flipud(f_hat(:));
    end
    
    f_hat_ts = stimDt*(-numStimFuture:numStimPast);
    
    % Only actually do the work of flipping S to the standard convention if
    % it's going to be outputted
    if nargout > 2 && fastestMethod < 3
        S = fliplr(S);
    end

end

function S = generateS_sparse(stim,sIdxs,numForward,numBack)

    S = zeros(length(sIdxs),numForward+numBack+1);

    for ii = 1:length(sIdxs)
        S(ii,:) = stim(sIdxs(ii)-numBack:sIdxs(ii)+numForward);
    end

end

function S = generateS_sparseTranspose(stim,sIdxs,numForward,numBack)

    S = zeros(numForward+numBack+1,length(sIdxs));

    for ii = 1:length(sIdxs)
        S(:,ii) = stim(sIdxs(ii)-numBack:sIdxs(ii)+numForward);
    end
    
    S = S';

end

% Note that this generates an S matrix with the columns flipped such that
% S_flipped == fliplr(S). This actually corresponds to the more standard
% convention of having responses predicting future stimuli on the left-hand
% side.
function S_flipped = generateS_dense(stim,sIdxs,numForward,numBack)

    row = zeros(1,numForward+numBack+1);
    row(1) = stim(1);
    toep = toeplitz(stim,row);

    S_flipped = toep(sIdxs+numForward,:);

end

function [ts] = predictSTime(numKernelElements,numResponses,numStimuli)
    tS  = 1.0e-8*numResponses*(60+numKernelElements);
    tST = 1.4e-8*numResponses*(20+numKernelElements);
    tD = (1.8e-8*numStimuli + 3e-9*numResponses)*numKernelElements;
    
    ts = [tS,tST,tD];
end