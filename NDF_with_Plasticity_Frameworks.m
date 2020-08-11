function NDF_with_Plasticity_Frameworks(a,lr,nTrialMax)
% clc;clear all;close all;    
datapath = ['../../data/FR_Curr_ring_RK4_distractor_with_Plasticity/' datestr(now,'yymmdd_HH_MM_') num2str(a) '_' num2str(lr) ];
mkdir(datapath)
%% load parameters
param = NDF_with_Plasticity_Parameters(a,lr,nTrialMax)
save([datapath,'/param.mat'],'-struct','param');

%% unpack Connectivity profile 
MEE = param.MEE;
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;

%% pack initial values
nx = param.N;
y0 = [0;              % Stimlus Current Strength
      0;              % Wipe Current Strength
      zeros(6*nx,1);  % 6*N state variables
      reshape(MEE,nx*nx,1) % E to E Connection Strength
      ]; 

%% load timings
TrialOn = param.TStimOn;
TrialEnd = param.TForgetOff;
dt_store = param.dt_store;
%% Solving ODE equations
options = odeset('RelTol',1e-3,'AbsTol',1e-5); 
disp(['Integration started at: ',datestr(now,'HH:MM:SS')])
nTB = param.nTrialBatch;
nTrial = nTB;
D = 1; % initially, set decoding error to maximum
DList = [];
while nTrial<=param.nTrialMax
    %% solve current batch
    [t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
        TrialOn(nTrial-nTB+1):dt_store:TrialEnd(nTrial),y0,options);
    %% unpack and save batch results
    nt = length(t);
    Mt = y(:,nx*6+3:end);MEEt = reshape(Mt,nt,nx,nx);
    t_readout = param.TDelayOff(1:nTB)-TrialOn(1);
    MEEt = MEEt(t_readout/dt_store,:,:);
    MEEt = permute(MEEt,[2 3 1]); % put time on 3rd dimention
    Rt = y(:,3:nx*6+2);Rt = reshape(Rt,nt,nx,6);
    RE = Rt(:,:,1)';RI = Rt(:,:,2)';SEE = Rt(:,:,3)';SIE = Rt(:,:,4)';SEI = Rt(:,:,5)';SII = Rt(:,:,6)'; % put time on 2nd dimention
    clear Rt;
    % Input_forget=y(:,2);
    % Input_stim = y(:,1); 
    % Input = Input_stim - Input_forget;
    save([datapath,'/results_' num2str(nTrial) '.mat'],'t','RE','RI','MEEt');
    %% plot and save
    h2=figure; %imagesc([RE RE1])
    tIndex = t>=TrialOn(nTrial-nTB+1) & t<TrialEnd(nTrial-nTB+10);
    subplot(2,1,1);imagesc(RE(:,tIndex),[0 50]);title(['trial ' num2str(nTrial-nTB+1) ' to ' num2str(nTrial-nTB+10) ])
    ylabel('position')
    tIndex = t>=TrialOn(nTrial-9) & t<TrialEnd(nTrial);
    subplot(2,1,2);imagesc(RE(:,tIndex),[0 50]);title(['trial ' num2str(nTrial-9) ' to ' num2str(nTrial)])
    ylabel('position')
    xlabel('Time')
    saveas(h2,[datapath,'/RE_' num2str(nTrial) '.jpg'])    
    %% calculate decoding performance
    theta = param.stimLoc_theta((nTrial-nTB+1):nTrial);
    t_readout = param.TDelayOff(1:nTB)-TrialOn(1);
    RE_readout = RE(:,t_readout/dt_store);
    nRE = poissrnd(RE_readout*0.2); % readout in a 200ms window. poisson spike count
    theta_hat = atan2(sin(param.x)*nRE,cos(param.x)*nRE);
    D = mean(1-cos(theta_hat-theta')); DList = [DList,D];
    %% update to next batch
    disp([num2str(nTrial) ' trials completed at: ',datestr(now,'HH:MM:SS'), '. D=',num2str(D)])
    y0 = y(end,:); % initial values
    nTrial = nTrial+nTB; 

end
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])

figure,plot(nTB:nTB:nTrial,DList)
xlabel('Trials')
ylabel('$\langle 1-\cos(\theta-\hat\theta) \rangle$','Interpreter','latex')


