function param = NDF_with_Plasticity_Parameters()
    %% Time constants
    % % neurons 
    param.TE = 20; % Excitatory population
    param.TI = 10; % Inhibitory population
    % % synapses
    param.TEE = 100; 
    param.TIE = 25;
    param.TEI = 10;
    param.TII = 10;
    % % stimulus filter
    param.Tinput = 100;
    
    %% Discretizing the space x
    nx = 64;
    dx = 2*pi/nx;
    x = -pi:dx:pi-dx; %periodic boundary 
    
    param.N = nx;
    param.dx= dx;
    param.x = x;
    
    %% Connectivity Profile
    J = 100;
    JEE = 1*J;JIE = 2*J;JEI = 1*J;JII = 2*J; 
    sigma_E = 0.2*pi; sigma_I = 0.1*pi; % wider excitatory synaptic projection
    f_E = 1*exp(-(x/sigma_E).^2);
    f_I = 1*exp(-(x/sigma_I).^2);


    MEE = zeros(nx,nx);
    MEI = zeros(nx,nx);

    for i = 1:nx
        MEE(i,:) = dx*circshift(f_E,[0 -pi/dx-1+i]);
        MEI(i,:) = dx*circshift(f_I,[0 -pi/dx-1+i]);
    end
    MIE = MEE;
    MII = MEI;

    MEE = JEE*MEE;
    MIE = JIE*MIE;
    MEI = JEI*MEI;
    MII = JII*MII;

    param.MEE = MEE;
    param.MIE = MIE;
    param.MEI = MEI;
    param.MII = MII;

    %% Transfer Function
    NE = 3;
    thE = 10;
    sigE = 100;
    maxfE = 100;
    qE = @(x) maxfE*(x-thE).^NE./(sigE^NE+(x-thE).^NE).*(x>thE);
    
    NI = 3;
    thI = 10;
    sigI = 100;
    maxfI = 100;
    qI = @(x) maxfI*(x-thI).^NI./(sigI^NI+(x-thI).^NI).*(x>thI);
    param.qE = qE;
    param.qI = qI;
%     param.qE = @(x) x.*(x>0);
%     param.qI = @(x) x.*(x>0);
%     param.qE = @(x) 7*((0.2*(x-1)+0.5).*(x<=1) + (2*(x-1)+0.5).*(x>1).*(x<=2) + (1*(x-2)+2.5).*(x>2).*(x<=14) + 14.5.*(x>14));
%     param.qI = @(x) x; 
    %% Perturbations
    a = .65;
    % %sharp local perturbation 
    % index_x = 0:dx:pi/8;
    % index = floor((index_x+pi)/dx)+1;
    % MEE0 = MEE; MEE(index,:) = a*MEE(index,:); 
    % param.MEE = MEE;
    % param.MEE_unperturbed = MEE0;
    % param.perturbation_type = 'local-rowwise(postsyn)';
    % param.perturbation_strength = a;    
%     MEE0 = MEE; MEE(:,index) = a*MEE(:,index); 
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'local-colwise(presyn)';
%     param.perturbation_strength = a;    
    % % smooth local perturbation
%     index_x = 0:dx:pi/8;
%     index = floor((index_x+pi)/dx)+1;
%     perturbation = 1 - (1-a)*exp(-((x-0.125*pi)/(pi/2)).^2);
%     param.perturbation = perturbation;
%     perturbation = repmat(perturbation',1,nx);
%     MEE0 = MEE; MEE = MEE.*perturbation;
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'local-rowwise(postsyn)';
%     param.perturbation_strength = a;    
%     param.perturbation_index = index_x;
% % global perturbation
    MEE0 = MEE; MEE = MEE*a;
    param.MEE = MEE;
    param.MEE_unperturbed = MEE0;
    param.perturbation_type = 'Global';
    param.perturbation_strength = a;
% % random perturbation
%     perturbation = 10.^(randn(nx)/1000);
%     MEE0 = MEE; MEE = MEE.*perturbation;
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'random';
%     previousResult = load("C:\Users\golde\Documents\Research\data\FR_Curr_ring_RK4_distractor_with_Plasticity\190917_19_47_\results.mat");
%     MEE0 = MEE;MEE = previousResult.MEEt(:,:,end);
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;

    %% External Input
    JEO = 2*J;
    IEO_init = 3.5*(exp(-(x/(pi/4)).^2)');
    
    param.IEo = JEO*IEO_init;
    param.IIo = 0;
    param.IEc = 15*ones(nx,1); 
    param.IIc = 0;
    
    %% simulation timing in milisecond

    T_on = 1000; % time between initialization and stimulus onset
    Tstim = 500; % stimulus presentation duration
    Tmemory = 4500; % delay duration
    
    nTrialBatch = 1e2; % number of trials after which to run a homeostasis training
    nBatch = 1e2; 
    nTrial=nTrialBatch*nBatch; % number of training trails

    THomeo = 1e6; % baseline duration to run homeostasis plasticity

    tTrial = T_on+Tstim+Tmemory;
    tBatch = tTrial*nTrialBatch + THomeo; 

    TrialOn_batch = 0:tTrial:tTrial*(nTrialBatch-1);
    
    TBatchOn = 0:tBatch:tBatch*(nBatch-1);

    TrialOn = repmat(TrialOn_batch,1,nBatch) + repelem(TBatchOn,nTrialBatch);
    TStimOn = TrialOn + T_on;
    TStimOff = TStimOn + Tstim;
    TrialOff = TStimOff + Tmemory;

    param.nTrial = nTrial; 
    param.nTrialBatch = nTrialBatch; 
    param.nBatch = nBatch;
    param.TrialOn = TrialOn;
    param.TrialOff = TrialOff;
    param.TStimOn   = TStimOn;
    param.TStimOff  = [TStimOff];
    param.TBatchOn = TBatchOn;
    param.THomeoOn = TBatchOn + tTrial*nTrialBatch;;
    param.TBatchOff = TBatchOn+tBatch;
    
    param.dt_store_delay = 50;
    param.dt_store_homeo = 1e3;

    %% randomize stimlus location
    stimLoc = randi(nx,nTrial,1)-nx/2; % in training period
    stimLoc_theta = stimLoc/nx*2*pi;
%     stimLoc_test = 0; 

    param.stimLoc = [stimLoc];
    param.stimLoc_theta = [ stimLoc_theta ];

    %% delay-period plasticity
    % x: nx by 1, x: post-syn, x': pre-syn
    param.fM_expr = '@(x,dx) ( (1) .* dx ) * x'' ';
    param.fM = eval(param.fM_expr);
    % % Plasticity 1/learning rate
    param.TJ_delay = 1e4;

    %% homeostatic plasticity
    param.TJ_homeostasis = 1e6;
    param.RE_target = 1; % determine from parameters above for now...
