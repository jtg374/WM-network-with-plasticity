function NDF_with_Plasticity_FillInMiddle(datapath,iTrialBegin,iTrialEnd)
    disp(datapath)

param = load([datapath,'param.mat']);
load([datapath,'results.mat']);
%% unpack Connectivity profile 
MEE = MEEt(:,:,iTrialBegin-1);
g = g_readout(:,iTrialBegin-1);
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;

%% pack initial values
nx = param.N;
np = param.np;
y0 = [0;              % Stimlus Current Strength
      0;              % Wipe Current Strength
      zeros(6*nx*np,1);  % 6*N state variables
      reshape(MEE,nx*nx,1); % E to E Connection Strength
      g % gains
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
for iTrial=iTrialBegin:iTrialEnd
    %% solve current batch
    [t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
        TrialOn(iTrial):dt_store:TrialEnd(iTrial),y0,options);
    %% unpack and save batch results
    nt = length(t);
    gt = y(:,nx*nx+nx*np*6+3:end)';g=gt(:,end);
    Mt = y(:,nx*np*6+3:nx*nx+nx*np*6+2);Mt = reshape(Mt,nt,nx,nx);
    MEE = Mt(end,:,:); MEE = squeeze(MEE);
    Rt = y(:,3:nx*np*6+2);Rt = reshape(Rt,nt,nx,np,6);Rt = permute(Rt,[2,3,1,4]);
    RE = Rt(:,:,:,1);RI = Rt(:,:,:,2);SEE = Rt(:,:,:,3);SIE = Rt(:,:,:,4);SEI = Rt(:,:,:,5);SII = Rt(:,:,:,6); 
    clear Rt;
    %% plot and save
    addpath('/gpfsnyu/home/jtg374/MATLAB/CubeHelix') 
    % if mod(iTrial,modTrial)==0 | ismember(iTrial,[1,2,5,10,20,50,100,200,500,1000,2000])
        save([datapath,'/FullData/results_' num2str(iTrial) '.mat'],'t','RE','RI','gt');
        h2=figure; %imagesc([RE RE1])
        imagesc(squeeze(RE(:,param.pNp(iTrial),:)),[0 50]);
        ylabel('neuron')
        xlabel('Time')
        colormap(cubehelix)
        saveas(h2,[datapath,'/ActFigures/RE_T_' num2str(iTrial) '.jpg'])
        saveas(h2,[datapath,'/ActFigures/RE_T_' num2str(iTrial) '.eps'],'epsc')
        h3=figure;
        imagesc(RE(:,:,end),[0 50])
        xlabel('stim position')
        ylabel('neuron')
        colormap(cubehelix)
        colorbar
        saveas(h3,[datapath,'/ActFigures/RE_X_' num2str(iTrial) '.jpg'])
        h4=figure;
        imagesc(diag(g)*MEE,[0 10])
        xlabel('pre-syn')
        xlabel('post-syn')
        colormap(cubehelix)
        colorbar
        saveas(h4,[datapath,'/ActFigures/g-MEE_' num2str(iTrial) '.jpg'])
    % end
    disp([num2str(iTrial) ' trials completed at: ',datestr(now,'HH:MM:SS')])

    %% update to next batch
    y0 = [  0;              % Stimlus Current Strength
            0;              % Wipe Current Strength
            zeros(6*nx*np,1);  % 6*N state variables
            reshape(MEE,nx*nx,1); % E to E Connection Strength
            g %gain
            ]; 
end
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])


