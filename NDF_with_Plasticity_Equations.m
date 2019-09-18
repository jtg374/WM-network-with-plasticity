function dy = NDF_with_Plasticity_Equations(t,y,param)

% unpack variables
N = param.N;
MEE = y(N*6+1:end);
MEE = reshape(MEE,N,N);MEE(MEE<0)=0;
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;
R = y(1:N*6);R = reshape(R,N,6);RE = R(:,1);RI = R(:,2);SEE = R(:,3);SIE = R(:,4);SEI = R(:,5);SII = R(:,6);

% load timing
iTrial    = sum(t>=param.TrialOn);
TStimOn   = param.TStimOn(iTrial);
TStimOff  = param.TStimOff(iTrial);

% additional parameters
IEc = param.IEc;
IIc = param.IIc;
% set stimulus location
shift=0;
if iTrial
    shift = param.stimLoc(iTrial);
end
IEo = circshift(param.IEo,shift);IIo = circshift(param.IIo,shift);

% determine stimulus strength
if t<TStimOn % before stimulus presentation, baseline, no input
    IStim = 0;
elseif t<TStimOff % during stimulus presentation
    IStim = 1 - exp(-(t-TStimOn)/param.Tinput);
else % delay period
    IStim = exp(-(t-TStimOff)/param.Tinput);
end

% transfer function
qE = param.qE;
qI = param.qI;

% main ode eqs
% % Neurons Populations and Synapses
dRe = 1./param.TE .*( -RE + qE(MEE*SEE - MEI*SEI + IEo*IStim + IEc));
dRi = 1./param.TI .*( -RI + qI(MIE*SIE - MII*SII + IIc));
dSee= 1./param.TEE.*(-SEE + RE);
dSie= 1./param.TIE.*(-SIE + RE);
dSei= 1./param.TEI.*(-SEI + RI);
dSii= 1./param.TII.*(-SII + RI);
% % Plasticity
K=1e3;dRe_ = dRe;dRe_(dRe>K)=K; % set an upper bound for plasticity
if t>TStimOff
    %
    fM = param.fM;
    dMEE= -1/param.TJ * fM(RE,dRe_);
    %
else 
    dMEE=zeros(N,N);
end

% pack variable derivatives
dy=[dRe;dRi;dSee;dSie;dSei;dSii;reshape(dMEE,N*N,1)];  