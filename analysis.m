datapath = 'C:\Users\golde\Documents\Research\data\FR_Curr_ring_RK4_distractor_with_Plasticity\200904_11_45_Global';
param = load([datapath,'\param.mat']);
load([datapath,'\results.mat']);
%% Decode
theta = param.stimLoc_theta;
t_readout = param.TDelayOff;
RE_readout = RE(:,:,t_readout/param.dt_store);
nRE = poissrnd(RE_readout*0.2); % readout in a 200ms window. poisson spike count
theta_hat = zeros(size(theta));
for ip = 1:param.np
    nREip = squeeze(nRE(:,ip,:));
    theta_hat(ip,:) = atan2(sin(param.x)*nREip,cos(param.x)*nREip);
end
D = 1-cos(theta_hat-theta);