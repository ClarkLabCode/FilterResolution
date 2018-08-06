%% script shows intuition for super-resolution

rng(2); % to get similar behavior each time

t = [0:1e5]; % 100 seconds in ms slots
tsamp = t(1:500);
f = exp(-tsamp/100).*sin(tsamp*2*pi/50);
f = f./sqrt(sum(f.^2)); % normalize L2

stimRate = 1/(3000); % 1 input every 3 seconds, on average
stim = double(rand(size(t))<stimRate); % binary stimulus variable
stimTimes = find(stim>0);

noiseamp = 0.003;
resp = filter(f,1,stim) + noiseamp*randn(size(t)); % noise + response

sampleInterval = 100; % sample every 100 ms
sampleTimes = [0:sampleInterval:max(t)];

% How far forward and back to calculate filter, in ms
extentForwards = 500;
extentBackwards = 100;

%% Zoomed out plot

figure(1);
clf;

subplot(3,1,1)
plot(t,stim);
xlim([stimTimes(1)-1e3 stimTimes(7)+1e3])

subplot(3,1,2); hold on;
for ii=1:length(sampleTimes)
    plot(sampleTimes(ii)*[1 1],[0 1],'k','linewidth',0.5);
end
xlim([stimTimes(1)-1e3 stimTimes(7)+1e3])

subplot(3,1,3);
plot(t,resp);
xlim([stimTimes(1)-1e3 stimTimes(7)+1e3])

for tt = stimTimes
    rectYExtent = 0.2;
    rectangle('Position',[tt-extentBackwards, -rectYExtent, extentBackwards+extentForwards, 2*rectYExtent]);
end

set(gcf, 'Renderer', 'painter');

%% Plot showing sampling around a single stimulus
figure(2);
clf;

subplot(3,1,1); hold on;
subplot(3,1,2); hold on;
subplot(3,1,3); hold on;
selectedStimTime = stimTimes(1);
snipStart = selectedStimTime-extentBackwards;
snipEnd = selectedStimTime+extentForwards;
plot([-extentBackwards:extentForwards],resp(snipStart:snipEnd));
for ii=1
    subplot(3,1,1);
    plot([-extentBackwards:extentForwards],stim(stimTimes(ii)+[-extentBackwards:extentForwards]));
    
    ff = find( (sampleTimes>=(stimTimes(ii)-extentBackwards)) .* (sampleTimes<=stimTimes(ii)+extentForwards) );
    subplot(3,1,2);
    for jj=1:length(ff)
        plot( (sampleTimes(ff(jj))-stimTimes(ii))*[1 1],[0 3]+(ii-1), 'k' );
    end
    
    subplot(3,1,3);
    scatter( sampleTimes(ff) - stimTimes(ii) , resp(sampleTimes(ff)+1) , 'k' ); 
    % bad coding: can use sampleTimes because they are also the index
    
end
subplot(3,1,2);
set(gca,'ylim',[0 3],'xlim',[-extentBackwards extentForwards]);

%% Plot with a few examples of stim and response with sampling
% f3 = figure(3);
% cmap = colormap('cool');
% cmapInd = round((size(cmap,1)-1)*sampleTimes/max(t)+1);
% 
% subplot(3,1,1)
% plot(t,stim);
% ylabel('Stimulus')
% 
% subplot(3,1,2); hold on;
% for ii=1:length(sampleTimes)
%     plot(sampleTimes(ii)*[1 1],[0 1],'color',cmap(cmapInd(ii),:));
% end
% ylabel('Sample')
% 
% subplot(3,1,3);
% plot(t,resp);
% ylabel('Response')
% 
% sampleTimesNearStim = [];
% extentForwards = 500;
% extentBackwards = 100;
% for stimNum = 1:length(stimTimes)
%     stimT = stimTimes(stimNum);
%     samplesNearThisStim = sampleTimes(((sampleTimes - stimT) < extentForwards) & ((stimT - sampleTimes) < extentBackwards));
%     sampleTimesNearStim = [sampleTimesNearStim samplesNearThisStim];
%     sampleTimesNearStimCell{stimNum} = samplesNearThisStim;
% end
% hold on;
% scatter(sampleTimesNearStim,resp(sampleTimesNearStim+1),'k');
% for ii = 1:length(sampleTimesNearStimCell)
%     plot(sampleTimesNearStimCell{ii},resp(sampleTimesNearStimCell{ii}+1),'k');
% end
% hold off;
% figure(f3);
% for ii=1:3
%     subplot(3,1,ii);
%     set(gca,'xlim',[7.5e3 1.25e4]);
% end


%% now make a plot showing how super-resolution works
figure(4);
clf;

subplot(3,1,1); hold on;
plot([-extentBackwards:extentForwards],stim(stimTimes(1)+[-extentBackwards:extentForwards]));

subplot(3,1,2); hold on;
subplot(3,1,3); hold on;
plot([-extentBackwards tsamp],[0 f],'LineWidth',1,'color',[0.5 .5 .5]);

numStimsToAverage = 20;
accumulatedSampleTimes = [];
accumulatedResps = [];
for ii=1:numStimsToAverage
    ff = find( (sampleTimes>=(stimTimes(ii)-extentBackwards)) .* (sampleTimes<=stimTimes(ii)+extentForwards) );
    subplot(3,1,2);
    for jj=1:length(ff)
        plot( (sampleTimes(ff(jj))-stimTimes(ii))*[1 1],[0.1 0.9]+(ii),'k');
        plot( (sampleTimes(ff(jj))-stimTimes(ii))*[1 1],[0.1 0.9],'k');
    end
   
    ax = gca;
    ax.YAxis.MinorTickValues = 1:numStimsToAverage;
    ax.YMinorGrid = 'on';
    
    subplot(3,1,3);
    scatter( sampleTimes(ff) - stimTimes(ii) , resp(sampleTimes(ff)+1),'k');
    
    accumulatedSampleTimes = [accumulatedSampleTimes sampleTimes(ff) - stimTimes(ii)];
    accumulatedResps = [accumulatedResps resp(sampleTimes(ff)+1)];
    % bad coding: can use sampleTimes because they are also the index    
end

% smoothed = fit(accumulatedSampleTimes',accumulatedResps','smoothingspline','SmoothingParam',0.05);
% plot(smoothed,'k')
% legend('hide')
% xlabel('time')
% ylabel('fit response')

subplot(3,1,2);
set(gca,'ylim',[0 numStimsToAverage+1]);
ylabel('stim number')

%% make new figures for cartoon of sampling

im = zeros(5,5);
im(3,4) = 1;
figure;
for ii=1:3
    subplot(2,3,ii); hold on;
    imagesc(im); colormap('gray');
    for jj = 1:5
        for kk = 1:5
            plot(jj,kk,'o','color',[.5 .5 .5]);
        end
        line([1 5],jj*[1 1],'linewidth',1,'color',[.5 .5 .5]);
        if jj<5
            line([5 1],[jj+1,jj],'linewidth',1,'color',[.5 .5 .5]);
        end
    end
    plot(4,3,'or');
    set(gca,'xtick',[],'ytick',[],'xlim',[0.5 5.5],'ylim',[0.5 5.5]);
    
end

subplot(2,1,2); hold on;
for ii=1:75
    if mod(ii,25) ~= 14
        plot(ii*[1 1],[0 0.75],'color',[0.5 .5 .5],'linewidth', 2);
    else
        plot(ii*[1 1],[0 1],'color',[0.8 .2 0.2],'linewidth',2);
    end
end
set(gca,'xtick',[],'ytick',[0 1],'xlim',[0 76]);

%% make neuron/neurite looking thing

[xi,yi] = meshgrid([1:100],[1:100]);
neur = 0.5*(1+tanh(5 - ((yi-50).^2 + (xi-70).^2)/20 ) );

figure;
subplot(2,3,1);
imagesc(neur); colormap('gray');
set(gca,'xtick',[],'ytick',[]);

%% okay, now do a Guassian input, and show the responses to get the cross-correlation

stim = randn(size(t));
resp = filter(f,1,stim);

figure; hold on;
subplot(3,1,1);
plot(stim(1000:2000));
set(gca,'xlim',[0 1000]);
subplot(3,1,2); hold on;
plot(resp(1000:2000));
plot([1:100:1001], resp(1000:100:2000),'ok');
for ii=1:100:1001
    plot(ii*[1 1],[3.5 4],'k-');
end
set(gca,'xlim',[0 1000]);

sampleTimesSel = sampleTimes( find( (sampleTimes > extentForwards+1) .* (sampleTimes < length(t)-extentBackwards - 1)));
respWeightedStim = zeros(length(sampleTimesSel),extentBackwards+extentForwards+1);
for ii=1:length(sampleTimesSel)
    % note, this is just an explicit computation of the cross-correlation
    temp = resp(sampleTimesSel(ii))*stim(sampleTimesSel(ii)+[-extentForwards:extentBackwards]);
    respWeightedStim(ii,:) = temp(end:-1:1);

end

% figure; hold on;
subplot(3,1,3); hold on;
% plot([-extentBackwards:extentForwards],respWeightedStim([1:100:end],:)','color',0.5*ones(1,3));
plot([-extentBackwards:extentForwards],mean(respWeightedStim,1),'k','linewidth',1);
plot([-100:0,1:length(f)],[zeros(1,101),f(:)'],'color',0.75*ones(3,1),'linewidth',2);
set(gca,'ylim',[-.2 0.25]);




