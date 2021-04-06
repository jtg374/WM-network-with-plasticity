function dy = NDF_with_Plasticity_Equations(t,y,param)

% unpack variables
N = param.N;
np = param.np;
g = y(N*N+N*np*6+3:end);
MEI = y(N*np*6+3:N*np*6+N*N+2);
MEI = reshape(MEI,N,N);
MEE = param.MEE;
MIE = param.MIE;
MII = param.MII;
R = y(3:N*np*6+2);R = reshape(R,N,np,6);RE = R(:,:,1);RI = R(:,:,2);SEE = R(:,:,3);SIE = R(:,:,4);SEI = R(:,:,5);SII = R(:,:,6);
IStim = y(1); IWipe=y(2);

% load timing
TStimOn   = param.TStimOn;
TStimOff  = param.TStimOff;
TDelayOff = param.TDelayOff;
TForgetOff= param.TForgetOff;
nTrial    = sum(t>=TStimOn);

% set stimulus location
shift=zeros(np,1);
if nTrial
    shift = param.stimLoc(:,nTrial);
    iS = param.pNp(nTrial);
end
IEo = zeros(N,np);IIo = zeros(N,np);
for ip=1:np
    IEo(:,ip) = circshift(param.IEo,shift(ip));
    IIo(:,ip) = circshift(param.IIo,shift(ip));
end
% transfer function
qE = param.qE;
qI = param.qI;

% main ode eqs
% % Neurons Populations and Synapses
dRe = 1./param.TE .*( -RE + qE(diag(g)*MEE.*(MEE>=0)*SEE - MEI*SEI + IEo*IStim)*(1-IWipe));
dRi = 1./param.TI .*( -RI + qI(MIE*SIE - MII*SII + IIo*IStim)*(1-IWipe));
dSee= 1./param.TEE.*(-SEE + RE);
dSie= 1./param.TIE.*(-SIE + RE);
dSei= 1./param.TEI.*(-SEI + RI);
dSii= 1./param.TII.*(-SII + RI);
% % External Stimilus
dIt = 1./param.Tinput .*( -IStim + sum(t>=TStimOn)   - sum(t>TStimOff)  );
dIw = 1./param.Tinput .*( -IWipe + sum(t>=(TDelayOff+0)) - sum(t>TForgetOff) );
% % Plasticity
% K=10/500;dRe_ = dRe;dRe_(dRe>K)=K; % set an upper bound for plasticity
if any( (t>TStimOff).* (t<TDelayOff) )
    %
    fM_IP = param.fM_IP;
    dMEI= fM_IP(RE(:,iS),RI(:,iS))* (1-IStim)-MEI.*(MEI<0);
    dg  = zeros(N,1);
    %
else 
    dMEI=-MEI.*(MEI<0);
    dg = zeros(N,1);
end

% pack variable derivatives
dy=[dIt;dIw;
    reshape(dRe,N*np,1);
    reshape(dRi,N*np,1);
    reshape(dSee,N*np,1);
    reshape(dSie,N*np,1);
    reshape(dSei,N*np,1);
    reshape(dSii,N*np,1);
    reshape(dMEI,N*N,1);
    dg];  