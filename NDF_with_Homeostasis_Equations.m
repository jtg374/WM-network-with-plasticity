function dy = NDF_with_Homeostasis_Equations(t,y,param)

% unpack variables
N = param.N;
MEE = y(N*6+1:end);
MEE = reshape(MEE,N,N);MEE(MEE<0)=0;
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;
R = y(1:N*6);R = reshape(R,N,6);RE = R(:,1);RI = R(:,2);SEE = R(:,3);SIE = R(:,4);SEI = R(:,5);SII = R(:,6);

% additional parameters
IEc = param.IEc;
IIc = param.IIc;


% transfer function
qE = param.qE;
qI = param.qI;

% main ode eqs
% % Neurons Populations and Synapses
dRe = 1./param.TE .*( -RE + qE(MEE*SEE - MEI*SEI + IEc));
dRi = 1./param.TI .*( -RI + qI(MIE*SIE - MII*SII + IIc));
dSee= 1./param.TEE.*(-SEE + RE);
dSie= 1./param.TIE.*(-SIE + RE);
dSei= 1./param.TEI.*(-SEI + RI);
dSii= 1./param.TII.*(-SII + RI);
% % Plasticity
TJ = param.TJ_homeostasis; % 1/learning rate
diff=param.RE_target - RE;
diff(diff>param.RE_target*0.8) = param.RE_target*0.8;
dMEE = 1/TJ * diag(diff) * MEE;

% pack variable derivatives
dy=[dRe;dRi;dSee;dSie;dSei;dSii;reshape(dMEE,N*N,1)];  
