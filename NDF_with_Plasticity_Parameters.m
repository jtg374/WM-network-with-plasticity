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
    % % Plasticity 1/learning rate
    param.TJ = 1e4;
    
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
    NE = 2;
    thE = 10;
    sigE = 40;
    maxfE = 100;
    qE = @(x) maxfE*(x-thE).^NE./(sigE^NE+(x-thE).^NE).*(x>thE);
    
    NI = 2;
    thI = 10;
    sigI = 40;
    maxfI = 100;
    qI = @(x) maxfI*(x-thI).^NI./(sigI^NI+(x-thI).^NI).*(x>thI);
%     param.qE = qE;
%     param.qI = qI;
    param.qE = @(x) x.*(x>0);
    param.qI = @(x) x.*(x>0);
    %% Perturbations
    a = .97;
    % %sharp local perturbation 
    index_x = 0:dx:pi/8;
    index = floor((index_x+pi)/dx)+1;
    MEE0 = MEE; MEE(index,:) = a*MEE(index,:); 
    param.MEE = MEE;
    param.MEE_unperturbed = MEE0;
    param.perturbation_type = 'local-rowwise(postsyn)';
    param.perturbation_strength = a;    
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
% % % global perturbation
%     MEE0 = MEE; MEE = MEE*a;
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'Global';
%     param.perturbation_strength = a;
% % random perturbation
%     perturbation = 10.^(randn(nx)/1000);
%     MEE0 = MEE; MEE = MEE.*perturbation;
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'random';
    % previousResult = load("C:\Users\golde\Documents\Research\data\FR_Curr_ring_RK4_distractor_with_Plasticity\190806_11_11_LinearLargePerturb\results.mat");
    % MEE0 = MEE;MEE = previousResult.MEEt(:,:,end);
    % param.MEE = MEE;
    % param.MEE_unperturbed = MEE0;

    %% External Input
    JEO = 2*J;
    IEO_init = 1.35*(exp(-(x/(pi/4)).^2)'+1*ones(nx,1));
    
    param.IEo = JEO*IEO_init;
    param.IIo = 0;
    
    %% simulation timing in milisecond

    T_on = 500; %ms
    Tstim = 500;
    Tmemory = 3000;
    Tforget = 1000;
    
    nTrial=1000; % number of trails
    tTrial = T_on+Tstim+Tmemory+Tforget; % length of a trial
    tMax = nTrial*tTrial;

    TStimOn = T_on:tTrial:tMax;

    param.nTrial = nTrial;
    param.TStimOn   = TStimOn;
    param.TStimOff  = TStimOn+Tstim;
    param.TDelayOff = TStimOn+Tstim+Tmemory;
    param.TForgetOff= TStimOn+Tstim+Tmemory+Tforget;
    param.Tmax = tMax;

    param.dt_store = 50;

    %% randomize stimlus location
    stimLoc = randi(nx,nTrial,1)-nx/2; % in training period
    stimLoc_theta = stimLoc/nx*2*pi;

    param.stimLoc = stimLoc;
    % param.stimLoc_theta = stimLoc_theta;

    %% additional parameters for plasticity
    % x: nx by 1, x: post-syn, x': pre-syn
    param.RE_target = 20;
    param.fM_expr = '@(x,dx,M,x_ref) -dx * x'' + 1e-3*diag(x_ref-x)*M ';
    param.fM = eval(param.fM_expr);
