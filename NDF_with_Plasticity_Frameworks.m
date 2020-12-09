function NDF_with_Plasticity_Frameworks(a,lrD,lrH,nTrial,r_target)
% clc;clear all;close all;    
datapath = ['/gpfsnyu/scratch/jtg374/WM_Plasticity/UniformPerturb/MultiplictiveHomeo/' 'UniformP' num2str(a*100) 'DLR' num2str(lrD) 'HLR' num2str(lrH) datestr(now,'_yymmdd_HH_MM') 'Trial' num2str(nTrial) ];
mkdir(datapath)
disp(datapath)
mkdir([datapath '/FullData'])
mkdir([datapath '/ActFigures'])
%% load parameters
param = NDF_with_Plasticity_Parameters(a,lrD,lrH,nTrial,r_target)
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
      ones(nx,1) % gains
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
RE_readout = zeros(nx,np,nTrial);
g_readout = zeros(nx,nTrial);
for iTrial=1:nTrial
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
    MEEt(:,:,iTrial) = MEE;    
    save([datapath,'/FullData/results_' num2str(iTrial) '.mat'],'t','RE','RI','gt');
    %% plot and save
    addpath('/gpfsnyu/home/jtg374/MATLAB/CubeHelix') 
    if mod(iTrial,100)==0 | ismember(iTrial,[1,2,5,10,20,50])
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
        saveas(h3,[datapath,'/ActFigures/RE_X_' num2str(iTrial) '.jpg'])
    end
    disp([num2str(iTrial) ' trials completed at: ',datestr(now,'HH:MM:SS')])
    RE_readout(:,:,iTrial) = RE(:,:,end);
    g_readout(:,iTrial) = gt(:,end);
    %% update to next batch
    y0 = [  0;              % Stimlus Current Strength
            0;              % Wipe Current Strength
            zeros(6*nx*np,1);  % 6*N state variables
            reshape(MEE,nx*nx,1) % E to E Connection Strength
            g]; 

end
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])

save([datapath,'/results.mat'],'RE_readout','MEEt','g_readout');
saveas(h3,[datapath,'/RE_X_' num2str(iTrial) '.jpg'])

h=figure;
plot(g_readout')
c = hsv(nx);
set(gca, 'ColorOrder',c, 'NextPlot','ReplaceChildren');
plot(g_readout')
xlabel('Trial')
ylabel('gain')

saveas(h,[datapath,'/gain.fig'])
saveas(h,[datapath,'/gain.jpg'])

