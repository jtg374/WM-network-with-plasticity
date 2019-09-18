function NDF_with_Plasticity_Framework(perturbation,nTrain)
%% load parameters
disp(perturbation)
param = NDF_with_Plasticity_Parameters(perturbation,nTrain)


%% unpack Connectivity profile 
MEE = param.MEE;
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;

%% output time resolution
dt_store = param.dt_store;

%% pack initial values
nx = param.N;
y0 = [0;              % Stimlus Current Strength
      0;              % Wipe Current Strength
      zeros(6*nx,1);  % 6*N state variables
      reshape(MEE,nx*nx,1) % E to E Connection Strength
      ]; 
%% load timing
TStimOn   = param.TStimOn;
TStimOff  = param.TStimOff;
TDelayOff = param.TDelayOff;
Tmax = param.Tmax;

%% Solving ODE equations
clear textprogressbar % will cause trouble if integrate with @odetpbar while textprogressbar is in environment
%options = odeset('RelTol',1e-3,'AbsTol',1e-5,'OutputFcn',@odetpbar); % will print progressbar
options = odeset('RelTol',1e-3,'AbsTol',1e-5);
disp(['Integration started at: ',datestr(now,'HH:MM:SS')])
[t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
    0:dt_store:Tmax,y0,options);
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])

nt = length(t);
Mt = y(:,nx*6+3:end);MEEt = reshape(Mt,nt,nx,nx);
MEEt = MEEt(TDelayOff/dt_store,:,:);
MEEt = permute(MEEt,[2 3 1]); % put time on 3rd dimention
Rt = y(:,3:nx*6+2);Rt = reshape(Rt,nt,nx,6);
RE = Rt(:,:,1)';RI = Rt(:,:,2)';SEE = Rt(:,:,3)';SIE = Rt(:,:,4)';SEI = Rt(:,:,5)';SII = Rt(:,:,6)'; % put time on 2nd dimention
clear Rt;
% Input_forget=y(:,2);
% Input_stim = y(:,1); 
% Input = Input_stim - Input_forget;
clear y

%% save results
datapath = ['../../Data/WM_Plasticity/XSRule',num2str(perturbation*100),'GlobalPerturb',datestr(now,'_yyyyddmm_HH_MM')];
mkdir(datapath)
save([datapath,'/param.mat'],'-struct','param');
save([datapath,'/results.mat'],'t','TDelayOff','RE','RI','MEEt');
MEE_final = MEEt(:,:,end);
save(['../../Data/WM_Plasticity/XSRule/FinalMEE',num2str(perturbation*100),'.mat'],'perturbation','MEE_final');
