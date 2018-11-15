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

%% Plot showing sampling around a single stimulus
figure(1);
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
    ylabel('stimulus')
    
    ff = find( (sampleTimes>=(stimTimes(ii)-extentBackwards)) .* (sampleTimes<=stimTimes(ii)+extentForwards) );
    subplot(3,1,2);
    for jj=1:length(ff)
        plot( (sampleTimes(ff(jj))-stimTimes(ii))*[1 1],[0 3]+(ii-1), 'k' );
    end
    
    subplot(3,1,3);
    scatter( sampleTimes(ff) - stimTimes(ii) , resp(sampleTimes(ff)+1) , 'k' );
    ylabel('response samples')
    xlabel('time after stimulus (ms)');
end
subplot(3,1,2);
set(gca,'ylim',[0 3],'xlim',[-extentBackwards extentForwards]);
ylabel('response');

%% Zoomed out plot

figure(2);
clf;

subplot(3,1,1)
plot(t,stim);
xlim([stimTimes(1)-1e3 stimTimes(7)+1e3])
ylabel('stimulus')

subplot(3,1,2); hold on;
for ii=1:length(sampleTimes)
    plot(sampleTimes(ii)*[1 1],[0 1],'k','linewidth',0.5);
end
xlim([stimTimes(1)-1e3 stimTimes(7)+1e3])
ylabel('response samples')

subplot(3,1,3);
plot(t,resp);
xlim([stimTimes(1)-1e3 stimTimes(7)+1e3])

for tt = stimTimes
    rectYExtent = 0.2;
    rectangle('Position',[tt-extentBackwards, -rectYExtent, extentBackwards+extentForwards, 2*rectYExtent]);
end
ylabel('response')
xlabel('time (ms)')

set(gcf, 'Renderer', 'painter');

%% now make a plot showing how super-resolution works
figure(3);
clf;

subplot(3,1,1); hold on;
plot([-extentBackwards:extentForwards],stim(stimTimes(1)+[-extentBackwards:extentForwards]));
ylabel('stimulus')

subplot(3,1,2); hold on;
ylabel('stimulus number')
title('response samples')
subplot(3,1,3); hold on;
plot([-extentBackwards tsamp],[0 f],'LineWidth',1,'color',[0.5 .5 .5]);
ylabel('response')
xlabel('time after stimulus (ms)')

numStimsToAverage = 20;
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
end

subplot(3,1,2);
set(gca,'ylim',[0 numStimsToAverage+1]);

%% okay, now do a Guassian input, and show the responses to get the cross-correlation
figure(4);
clf;
hold on;

stim = randn(size(t));
resp = filter(f,1,stim);

subplot(3,1,1);
plot(stim(1000:2000));
set(gca,'xlim',[0 1000]);
ylabel('stimulus')

subplot(3,1,2); hold on;
plot(resp(1000:2000));
plot([1:100:1001], resp(1000:100:2000),'ok');
for ii=1:100:1001
    plot(ii*[1 1],[3.5 4],'k-');
end
set(gca,'xlim',[0 1000]);
ylabel('response samples')

sampleTimesSel = sampleTimes( find( (sampleTimes > extentForwards+1) .* (sampleTimes < length(t)-extentBackwards - 1)));
respWeightedStim = zeros(length(sampleTimesSel),extentBackwards+extentForwards+1);
for ii=1:length(sampleTimesSel)
    % note, this is just an explicit computation of the cross-correlation
    temp = resp(sampleTimesSel(ii))*stim(sampleTimesSel(ii)+[-extentForwards:extentBackwards]);
    respWeightedStim(ii,:) = temp(end:-1:1);

end

subplot(3,1,3); hold on;
plot([-extentBackwards:extentForwards],mean(respWeightedStim,1),'k','linewidth',1);
plot([-100:0,1:length(f)],[zeros(1,101),f(:)'],'color',0.75*ones(3,1),'linewidth',2);
set(gca,'ylim',[-.2 0.25]);
ylabel('filter weight')
xlabel('time before response (ms)')