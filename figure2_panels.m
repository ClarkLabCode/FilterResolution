%% script makes a short tutorial to show how to do cross-correlations correctly

rng(2); % do same thing every time

t = [0:100];
t1 = [-0.5:100.5];
stim = randn(size(t));
t2 = reshape([t1;t1],[1 2*length(t1)]);
t2 = t2(2:end-1);
for ii=1:length(t)
    stim2(1+2*(ii-1)+[0 1]) = stim(ii);
end
ti = [0:30];
B = [0, exp(-ti/10)/10];
A = 1;
resp = filter(B,A,stim);

sampint = 10;
sampi = find(mod(t,sampint)==0);

respnan = resp;
respnan(find(mod(t,sampint)))=nan;

respinterp = interpolateNan(respnan);

%% plot set up, with discrete vectors of stimulus and response

figure;
subplot(3,1,1); hold on;

surface([t2;t2],[stim2;stim2],zeros(2,length(t2)),[stim2;stim2],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);

plot([-0.5 100.5],[0 0],'k-');
set(gca,'xlim',20+[-0.5 50.5],'ylim',4*[-1 1]);


subplot(3,1,2); hold on;

surface([t;t],[resp;resp],zeros(2,length(resp)),[resp;resp],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);


plot(t(sampi),resp(sampi),'-ok');
for ii=1:length(sampi)
    line([1 1]*t(sampi(ii)),[0.5 0.6],'color',[0 0 0],'linew',2);
end
set(gca,'xlim',20+[-0.5 50.5],'ylim',0.6*[-1 1]);


subplot(9,1,7);
pcolor([stim(21:71) 1; stim(20:70) 1]);
colormap(copper);
set(gca,'dataa',[1 1 1],'xtick',[],'ytick',[]);
set(gca,'clim',[min(stim) max(stim)]);
box off;

subplot(9,1,8);
pcolor([respnan(21:71) 1; respnan(20:70) 1]);
colormap(copper);
set(gca,'dataa',[1 1 1],'xtick',[],'ytick',[]);
set(gca,'clim',[min(resp) max(resp)]);
box off;

subplot(9,1,9);
pcolor([respinterp(21:71) 1; respinterp(20:70) 1]);
colormap(copper);
set(gca,'dataa',[1 1 1],'xtick',[],'ytick',[]);
set(gca,'clim',[min(resp) max(resp)]);
box off;


%% plot snippets of stimulus and response

figure;
subplot(5,1,[1 3]); hold on;

% make the response have the same scale as stim, since colormap will match
resp2 = (resp-min(resp))/(max(resp)-min(resp))*(max(stim)-min(stim)) + min(stim);
xloc = 1:42;

for ii=3:8
    sloc = stim([-20:20]+sampi(ii));
    pcolor([xloc;xloc],[zeros(1,42);ones(1,42)]+ii*2,[sloc,1;sloc,1]);
    pcolor([45 46; 45 46],[0 0; 1 1] + ii*2,ones(2)*resp2(sampi(ii)));
end
colormap(copper);
set(gca,'dataa',[1 1 1],'xtick',[],'ytick',[]);
set(gca,'clim',[min(resp2) max(resp2)]);
box off;


stimlong = randn(1,1e4);
resplong = filter(B,A,stimlong);
samplongi = [21:10:length(stimlong)-20];

% create an interpolated version of the sampled resplong
resplongInterpolated = nan(size(resplong));
resplongInterpolated(samplongi)=resplong(samplongi);
resplongInterpolated = interpolateNan(resplongInterpolated);

for ii=1:length(samplongi)
    xc(ii,:) = resplong(samplongi(ii))*stimlong(samplongi(ii)+[-20:20]);
end
xcm = mean(xc,1);

for ii=100:length(stimlong)-100
    xc2(ii,:) = resplongInterpolated(ii)*stimlong(ii+[-20:20]);
end
xcm2 = mean(xc2,1);

subplot(5,1,4); hold on;
pcolor([xloc;xloc],[zeros(1,42);ones(1,42)],[xcm,1;xcm,1]);
pcolor([45 46; 45 46],[0 0; 1 1],zeros(2));
colormap(copper);
set(gca,'dataa',[1 1 1],'xtick',[],'ytick',[]);
set(gca,'clim',[min(xcm) max(xcm)]);
box off;

subplot(5,1,5);
plot([-20:20],xcm,'.k-'); hold on;
plot([-20:20],xcm2,'.m-');

% add groundtruth curve
groundTruth = filter(B,A,[zeros(1,20),1,zeros(1,20)]);
groundTruth = groundTruth(end:-1:1);
plot([-20:20],groundTruth,'-','color',[0.5 0.5 0.5],'linew',2);

xlabel('time relative to response');
ylabel('x-c with stim');
set(gca,'xlim',[-21,20+9]);



