function NDF_with_Plasticity_Frameworks(a,lrD,lrH,nTrialMax)
% clc;clear all;close all;    
datapath = ['/gpfsnyu/scratch/jtg374/WM_Plasticity_parameterSearch/' datestr(now,'yymmdd_HH_MM_') num2str(a) '_' num2str(lrD) ];
mkdir(datapath)
%% load parameters
param = NDF_with_Plasticity_Parameters(a,lrD,lrH,nTrialMax)
save([datapath,'/param.mat'],'-struct','param');

%% unpack Connectivity profile 
MEE = param.MEE;
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
TrialEnd = param.TDelayOff;
dt_store = param.dt_store;
%% Solving ODE equations
options = odeset('RelTol',1e-3,'AbsTol',1e-5); 
disp(['Integration started at: ',datestr(now,'HH:MM:SS')])
nTrial = param.nTrial;
D = 1; % initially, set decoding error to maximum
DList = [];
MEEt = zeros(nx,nx,nTrial);
RE_readout = zeros(nx,nTrial);
for iTrial=1:nTrial
    %% solve current batch
    [t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
        TrialOn(iTrial):dt_store:TrialEnd(iTrial),y0,options);
    %% unpack and save batch results
    nt = length(t);
    Mt = y(:,nx*np*6+3:end);Mt = reshape(Mt,nt,nx,nx);
    MEE = Mt(end,:,:); MEE = squeeze(MEE);
    Rt = y(:,3:nx*np*6+2);Rt = reshape(Rt,nt,nx,np,6);Rt = permute(Rt,[2,3,1,4]);
    RE = Rt(:,:,:,1);RI = Rt(:,:,:,2);SEE = Rt(:,:,:,3);SIE = Rt(:,:,:,4);SEI = Rt(:,:,:,5);SII = Rt(:,:,:,6); 
    clear Rt;
    % Input_forget=y(:,2);
    % Input_stim = y(:,1); 
    % Input = Input_stim - Input_forget;
    save([datapath,'/results_' num2str(iTrial) '.mat'],'t','RE','RI');
    %% apply homeostatic rule
    r_target = 15;
    r_mean = mean(mean(RE,3),2);
    MEE = MEE + diag(r_target-r_mean)*MEE*tTrial*lrH;
    MEEt(:,:,iTrial) = MEE;
    %% plot and save
    if mod(iTrial,100)==0
        h2=figure; %imagesc([RE RE1])
        imagesc(squeeze(RE(:,param.pNp(iTrial),:)),[0 50]);
        ylabel('position')
        xlabel('Time')
        saveas(h2,[datapath,'/RE_' num2str(iTrial) '.jpg'])    
    end
    disp([num2str(iTrial) ' trials completed at: ',datestr(now,'HH:MM:SS'), '. R_bar=',num2str(mean(r_mean))])
    RE_readout(:,iTrial) = RE(:,end);
    %% update to next batch
    y0 = [  0;              % Stimlus Current Strength
            0;              % Wipe Current Strength
            zeros(6*nx*np,1);  % 6*N state variables
            reshape(MEE,nx*nx,1) % E to E Connection Strength
            ]; 
    nTrial = nTrial+1; 

end
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])
