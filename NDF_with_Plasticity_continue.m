function NDF_with_Plasticity_continue(datapath,nTrialAdd)
    disp(datapath)

param = load([datapath,'param.mat']);
nx = param.N;
np = param.np;

T_on = param.TStimOn(1);
Tstim = param.TStimOff(1)-param.TStimOn(1);
Tmemory = param.TDelayOff(1) - param.TStimOff(1);
Tforget = param.TForgetOff(1) - param.TDelayOff(1);
tTrial = T_on+Tstim+Tmemory+Tforget; % length of a trial
tMin = param.Tmax;
tMax = param.Tmax + nTrialAdd*tTrial;
param.Tmax = tMax;

TStimOn = T_on:tTrial:tMax;

nTrialOld = param.nTrial;
param.nTrial = param.nTrial + nTrialAdd
param.TStimOn   = TStimOn;
param.TrialOn   = TStimOn-T_on;
param.TStimOff  = TStimOn+Tstim;
param.TDelayOff = TStimOn+Tstim+Tmemory;
param.TForgetOff= TStimOn+Tstim+Tmemory+Tforget;

stimLoc = randi(floor(nx/np),np,nTrialAdd); % random location in each group (1-4)
stimLoc = stimLoc + (0:floor(nx/np):(nx-1))'; % add level to each group
stimLoc = stimLoc - nx/2; % center to zero
stimLoc_theta = stimLoc/nx*2*pi;
pNp = randi(np,nTrialAdd,1);

param.stimLoc = [param.stimLoc stimLoc];
param.stimLoc_theta = [param.stimLoc_theta stimLoc_theta];
param.pNp = [param.pNp(:,1);pNp];

load([datapath,'results.mat']);
%% unpack Connectivity profile 
MEE = MEEt(:,:,end);
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;

%% pack initial values
nx = param.N;
np = param.np;
y0 = [0;              % Stimlus Current Strength
      0;              % Wipe Current Strength
      zeros(6*nx*np,1);  % 6*N state variables
      reshape(MEE,nx*nx,1) % E to E Connection Strength
      ]; 

%% load timings
tTrial = param.tTrial;
TrialOn = param.TrialOn;
T_on = param.TStimOn(1);
TrialEnd = param.TDelayOff;
dt_store = param.dt_store;
%% load plasticity
lrD = param.LearningRateDifferential;    
lrH = param.LearningRateHomeostatic;
r_target = param.r_target;


%% Solving ODE equations
options = odeset('RelTol',1e-3,'AbsTol',1e-5); 
disp(['Integration started at: ',datestr(now,'HH:MM:SS')])
new.MEEt = zeros(nx,nx,nTrialAdd);
new.RE_readout = zeros(nx,np,nTrialAdd);
MEEt = cat(3,MEEt,new.MEEt);
RE_readout = cat(3,RE_readout,new.RE_readout);
clear new
for iTrial=(nTrialOld+1):param.nTrial
    %% solve current batch
    [t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
        TrialOn(iTrial):dt_store:TrialEnd(iTrial),y0,options);
    %% unpack and save batch results
    nt = length(t);
    Rt = y(:,3:nx*np*6+2);Rt = reshape(Rt,nt,nx,np,6);Rt = permute(Rt,[2,3,1,4]);
    RE = Rt(:,:,:,1);RI = Rt(:,:,:,2);SEE = Rt(:,:,:,3);SIE = Rt(:,:,:,4);SEI = Rt(:,:,:,5);SII = Rt(:,:,:,6); 
    clear Rt;
    % Input_forget=y(:,2);
    % Input_stim = y(:,1); 
    % Input = Input_stim - Input_forget;
    %% apply homeostatic rule
    r_target = param.r_target;
    r_mean = mean(mean(RE,3),2);
    MEE = MEE + diag(r_target-r_mean)*MEE*tTrial*lrH;
    MEEt(:,:,iTrial) = MEE;
    %% apply XS rule
    dr = diff(RE,1,3);
    K = 10/500*dt_store; dr(dr>K) = K;
    for it = (T_on/dt_store):(nt-1)
        MEE = MEE - lrD*dr(:,:,it)*RE(:,:,it)'/nx;
    end
    MEEt(:,:,iTrial) = MEE;    
    %% plot and save
    addpath('/gpfsnyu/home/jtg374/MATLAB/CubeHelix') 
    if mod(iTrial,100)==0 | ismember(iTrial,[1,2,5,10,20,50])
        save([datapath,'/FullData/results_' num2str(iTrial) '.mat'],'t','RE','RI');
        h2=figure; %imagesc([RE RE1])
        imagesc(squeeze(RE(:,param.pNp(iTrial),:)),[0 50]);
        ylabel('neuron')
        xlabel('Time')
        colormap(cubehelix)
        saveas(h2,[datapath,'/ActFigures/RE_T_' num2str(iTrial) '.jpg'])
        h3=figure;
        imagesc(RE(:,:,end),[0 50])
        xlabel('stim position')
        ylabel('neuron')
        colormap(cubehelix)
        colorbar
        saveas(h3,[datapath,'/ActFigures/RE_X_' num2str(iTrial) '.jpg'])
    end
    disp([num2str(iTrial) ' trials completed at: ',datestr(now,'HH:MM:SS'), '. R_bar=',num2str(mean(r_mean))])
    RE_readout(:,:,iTrial) = RE(:,:,end);
    %% update to next batch
    y0 = [  0;              % Stimlus Current Strength
            0;              % Wipe Current Strength
            zeros(6*nx*np,1);  % 6*N state variables
            reshape(MEE,nx*nx,1) % E to E Connection Strength
            ]; 
end
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])

%% save results
disp(datapath)
save([datapath,'/results.mat'],'RE_readout','MEEt','-v7.3');
save([datapath,'/param.mat'],'-struct','param');
saveas(h3,[datapath,'/RE_X_' num2str(iTrial) '.jpg'])