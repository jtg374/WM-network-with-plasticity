%% Load Data 
datapath = '/gpfsnyu/scratch/jtg374/WM_Plasticity/randomPerturb_N256_200409_09_38/';
data1 = load([datapath,'results.mat']);
data1.param = load([datapath,'param.mat'])
%%
N = data1.param.N;
Nt = length(data1.param.TStimOn);
FFt = zeros(N,N,Nt);
FFEEt = zeros(N,N,Nt);
for ii = 1:Nt
    it = ii;
    FFt(:,:,ii) = fft2((data1.MEEt(:,:,it)-eye(N)) - data1.param.MEI/ (data1.param.MII + eye(N)) * data1.param.MIE);
end
save([datapath,'FF.mat'],'FFt')
data1.FF = load([datapath,'FF.mat']);

%% Page settings, position, font, etc
% create a A4 figure
close all
fig = figure();
fig.Units = 'points';
fig.Position = [0,0,595,842];
% create subplots
fSize = 10; % default font size

width = 0.15;
height = 0.09*0.8;
ax11 = axes();
ax11.Position = [0.4+width*0.1,.9-height,width*0.8,height]; %left,bottom,width,height
ax11.ActivePositionProperty = 'position';
ax11.FontSize = fSize;
ax12 = axes();
ax12.Position = [0.4+width*1.1,.9-height,width*0.8,height]; %left,bottom,width,height
ax12.ActivePositionProperty = 'position';
ax12.FontSize = fSize;
ax13 = axes();
ax13.Position = [0.4+width*2.1,.9-height,width*0.8,height]; %left,bottom,width,height
ax13.ActivePositionProperty = 'position';
ax13.FontSize = fSize;
cbPosition = [0.4+width*2.95,.9-height,0.01,height];


left = 0.4;
top = .9-height/0.8-0.04;
height = width*0.8;
ax21 = axes();
ax21.Position = [left+width*0.1,top-height,width*0.8,height];
ax21.ActivePositionProperty = 'position';
ax21.FontSize = fSize;
ax22 = axes();
ax22.Position = [left+width*1.1,top-height,width*0.8,height];
ax22.ActivePositionProperty = 'position';
ax22.FontSize = fSize;
ax23 = axes();
ax23.Position = [left+width*2.1,top-height,width*0.8,height];
ax23.ActivePositionProperty = 'position';
ax23.FontSize = fSize;
lgPosition2 = [left+width*2.95,top-height,width*0.8,0.06];

top = top-height/0.8-0.04;
height = 0.12*0.8;
ax31 = axes();
ax31.Position = [left+width*0.1,top-height,width*2.8,height];
ax31.ActivePositionProperty = 'position';
ax31.FontSize = fSize;
ax31.ColorOrder=summer(3);


%% Analysis and Plot
% nTrial = data1.param.nTrialTrain
N = data1.param.N;
x = data1.param.x;
nTrialPlot = 5;
nT = data1.param.nTrial;
tTrial = data1.param.TStimOn(2)-data1.param.TStimOn(1);
ntTrial = tTrial/data1.param.dt_store;
amp_range = [0 30];
f_m = 4;

iTrial =1;
axes(ax11)
RE = data1.RE(:,data1.t>=data1.param.TStimOn(iTrial)&data1.t<data1.param.TDelayOff(iTrial+nTrialPlot-1));
t = data1.t(data1.t>=data1.param.TStimOn(iTrial)&data1.t<data1.param.TDelayOff(iTrial+nTrialPlot-1))/tTrial;
imagesc(RE,[0,50])
xticks(1:ntTrial*2:(ntTrial*nTrialPlot))
xticklabels((1:2:nTrialPlot)+iTrial-1)
yticks([0,N/2,N])
yticklabels({'-\pi','0','\pi'})
ylabel('position')
box off

axes(ax21)
FF = fftshift(abs(data1.FF.FFt(:,:,iTrial)));
imagesc(FF((N/2+1-f_m):(N/2+1+f_m),(N/2+1-f_m):(N/2+1+f_m)),amp_range)
xticks(1:f_m/2:f_m*2+1)
yticks(1:f_m/2:f_m*2+1)
xticklabels(-f_m:f_m/2:f_m)
yticklabels(-f_m:f_m/2:f_m)
xlabel('\omega_j')
ylabel('\omega_i')
gca.FontSize = fSize;

iTrial =50;
axes(ax12)
box off
RE = data1.RE(:,data1.t>=data1.param.TStimOn(iTrial)&data1.t<data1.param.TDelayOff(iTrial+nTrialPlot-1));
t = data1.t(data1.t>=data1.param.TStimOn(iTrial)&data1.t<data1.param.TDelayOff(iTrial+nTrialPlot-1))/tTrial;
imagesc(RE,[0,50])
xticks(1:ntTrial*2:(ntTrial*nTrialPlot))
xticklabels((1:2:nTrialPlot)+iTrial-1)
yticks([])

axes(ax22)
FF = fftshift(abs(data1.FF.FFt(:,:,iTrial)));
imagesc(FF((N/2+1-f_m):(N/2+1+f_m),(N/2+1-f_m):(N/2+1+f_m)),amp_range)
xticks(1:f_m/2:f_m*2+1)
yticks(1:f_m/2:f_m*2+1)
xticklabels(-f_m:f_m/2:f_m)
yticklabels(-f_m:f_m/2:f_m)
xlabel('\omega_j')
gca.FontSize = fSize;


iTrial =450;
axes(ax13)
box off
RE = data1.RE(:,data1.t>=data1.param.TStimOn(iTrial)&data1.t<data1.param.TDelayOff(iTrial+nTrialPlot-1));
t = data1.t(data1.t>=data1.param.TStimOn(iTrial)&data1.t<data1.param.TDelayOff(iTrial+nTrialPlot-1))/tTrial;
imagesc(RE,[0,50])
xticks(1:ntTrial*2:(ntTrial*nTrialPlot))
xticklabels((1:2:nTrialPlot)+iTrial-1)
yticks([])
cb = colorbar();
cb.Position = cbPosition;
cb.Label.String = 'activity';
xlabel('Trial')


axes(ax23)
FF = fftshift(abs(data1.FF.FFt(:,:,iTrial)));
imagesc(FF((N/2+1-f_m):(N/2+1+f_m),(N/2+1-f_m):(N/2+1+f_m)),amp_range)
xticks(1:f_m/2:f_m*2+1)
yticks(1:f_m/2:f_m*2+1)
xticklabels(-f_m:f_m/2:f_m)
yticklabels(-f_m:f_m/2:f_m)
xlabel('\omega_j')
gca.FontSize = fSize;

axes(ax31)
offDiagVar = zeros(nT,1);
for it = 1:nT
    FF = data1.FF.FFt(:,:,it);
    FF = fftshift(abs(FF));
    FF = FF((N/2+1-f_m):(N/2+1+f_m),(N/2+1-f_m):(N/2+1+f_m));
    offDiagVar(it) = mean(abs(FF(~fliplr(eye(f_m*2+1)))).^2);
end
plot(offDiagVar)
xlabel('Trial')
ylabel('off-anti-diag variance')
ylim([0 max(offDiagVar)])
box off
%%
saveas(fig,[datapath,'fig5.pdf'])
saveas(fig,[datapath,'fig5.eps'])