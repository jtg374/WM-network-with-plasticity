function RE=NDF_with_Plasticity_1TrialTest(a,b)
% clc;clear all;close all;    
%% load parameters
param = NDF_with_Plasticity_Parameters(a,b,0,0,1,15)
% save([datapath,'/param.mat'],'-struct','param');

%% unpack Connectivity profile 
MEE = param.MEE;
MEI = param.MEI;
MIE = param.MIE;

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
iTrial=1;
    [t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
        TrialOn(iTrial):dt_store:TrialEnd(iTrial),y0,options);
    nt = length(t);
    Rt = y(:,3:nx*np*6+2);Rt = reshape(Rt,nt,nx,np,6);Rt = permute(Rt,[2,3,1,4]);
    RE = Rt(:,:,:,1);RI = Rt(:,:,:,2);SEE = Rt(:,:,:,3);SIE = Rt(:,:,:,4);SEI = Rt(:,:,:,5);SII = Rt(:,:,:,6); 
    clear Rt;
%%
x = param.x;
figure,imshow(squeeze(RE(:,randi(nx),:)),[],'InitialMagnification','fit')
figure;hold on; plot(x,RE(:,nx/2,end)); plot(x,RE(:,nx/2,10))
