function NDF_with_Plasticity_continue(datapath,nTrial)

param = load([datapath,'param.mat']);
nx = param.N;
np = param.np;

T_on = param.TStimOn(1);
Tstim = param.TStimOff(1)-param.TStimOn(1);
Tmemory = param.TDelayOff(1) - param.TStimOff(1);
Tforget = param.TForgetOff(1) - param.TDelayOff(1);
tTrial = T_on+Tstim+Tmemory+Tforget; % length of a trial
tMin = param.Tmax;
tMax = param.Tmax + nTrial*tTrial;
param.Tmax = tMax;

TStimOn = T_on:tTrial:tMax;

nTrialOld = param.nTrial;
param.nTrial = param.nTrial + nTrial
param.TStimOn   = TStimOn;
param.TStimOff  = TStimOn+Tstim;
param.TDelayOff = TStimOn+Tstim+Tmemory;
param.TForgetOff= TStimOn+Tstim+Tmemory+Tforget;

stimLoc = randi(floor(nx/np),np,nTrial); % random location in each group (1-4)
stimLoc = stimLoc + (0:floor(nx/np):(nx-1))'; % add level to each group
stimLoc = stimLoc - nx/2; % center to zero
stimLoc_theta = stimLoc/nx*2*pi;
pNp = randi(np,nTrial,1);

param.stimLoc = [param.stimLoc stimLoc];
param.stimLoc_theta = [param.stimLoc_theta stimLoc_theta];
param.pNp = [param.pNp(:,1);pNp];

load([datapath,'results.mat']);
%% unpack Connectivity profile 
MEE = MEEt(:,:,end);
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;

%% output time resolution
dt_store = param.dt_store;

%% pack initial values
nx = param.N;
np = param.np;
y0 = [0;              % Stimlus Current Strength
      0;              % Wipe Current Strength
      zeros(6*nx*np,1);  % 6*N state variables
      reshape(MEE,nx*nx,1) % E to E Connection Strength
      ]; 
%% load timing
TStimOn   = param.TStimOn;
TStimOff  = param.TStimOff;
TDelayOff = param.TDelayOff;
Tmax = param.Tmax;

%% Solving ODE equations
clear textprogressbar % will cause trouble if integrate with @odetpbar while textprogressbar is in environment
options = odeset('RelTol',1e-3,'AbsTol',1e-5,'OutputFcn',@odetpbar); % will print progressbar
disp(['Integration started at: ',datestr(now,'HH:MM:SS')])
[t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
    tMin:dt_store:Tmax,y0,options);
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])
%%
nt = length(t);
Mt = y(:,nx*np*6+3:end);new.MEEt = reshape(Mt,nt,nx,nx);
new.MEEt = new.MEEt(((T_on+Tstim+Tmemory):tTrial:(nTrial*tTrial))/dt_store,:,:);
new.MEEt = permute(new.MEEt,[2 3 1]); % put time on 3rd dimention
Rt = y(:,3:nx*np*6+2);Rt = reshape(Rt,nt,nx,np,6);Rt = permute(Rt,[2,3,1,4]);
RE = Rt(:,:,:,1);RI = Rt(:,:,:,2);SEE = Rt(:,:,:,3);SIE = Rt(:,:,:,4);SEI = Rt(:,:,:,5);SII = Rt(:,:,:,6); 

%% Figures
close all

%% 

h2=figure(2); %imagesc([RE RE1])
tIndex = t>=TStimOn(nTrialOld) & t<TDelayOff(nTrialOld+9);
subplot(2,1,1);imagesc(squeeze(RE(:,1,tIndex)),[0 50]);title('first 10 trials')
ylabel('position (80\theta / 2\pi)','FontSize',10)
tIndex = t>=param.TStimOn(end-9) & t<TDelayOff(end);
subplot(2,1,2);imagesc(squeeze(RE(:,1,tIndex)),[0 50]);title('last 10 trials')
ylabel('position (80\theta / 2\pi)','FontSize',10)
xlabel('Time (a.u.)','FontSize',14)
saveas(h2,[datapath,'/2-' num2str(param.nTrial) '.fig'])
saveas(h2,[datapath,'/2-' num2str(param.nTrial) '.jpg'])

%% save results
disp(datapath)
save([datapath,'/resultsFull' num2str(param.nTrial) '.mat'],'t','RE','RI','-v7.3');
new.RE_readout = RE(:,:,((T_on+Tstim+Tmemory):tTrial:(nTrial*tTrial))/dt_store);
MEEt = cat(3,MEEt,new.MEEt);
RE_readout = cat(3,RE_readout,new.RE_readout);
save([datapath,'/results.mat'],'t','TDelayOff','RE_readout','MEEt','-v7.3');
save([datapath,'/param.mat'],'-struct','param');