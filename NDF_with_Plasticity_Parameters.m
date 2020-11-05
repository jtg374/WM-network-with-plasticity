function param = NDF_with_Plasticity_Parameters(a,lrD,lrH,nTrialMax,r_target)
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
    np = nx;
    
    param.N = nx;
    param.dx= dx;
    param.x = x;
    param.np = np;
    
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
    param.qE = @(x) x.*(x>0) - (x-maxfE).*(x>maxfE);
    param.qI = @(x) x.*(x>0) - (x-maxfI).*(x>maxfI);
    %% Perturbations
%     a = 0.9;
    % %sharp local perturbation 
%     index_x = 0:dx:pi/8;
%     index = floor((index_x+pi)/dx)+1;
%     MEE0 = MEE; MEE(index,:) = a*MEE(index,:); 
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'local-rowwise(postsyn)';
%     param.perturbation_strength = a;    
% % combine global perturb
%     MEE = MEE*.9;
%     param.MEE = MEE;
%     param.perturbation_type = 'local-rowwise(postsyn)+global';
%     param.perturbation_strength_global = a*.9;    
    % %sharp local perturbation 
    % index_x = 0:dx:pi/8;
    % index = floor((index_x+pi)/dx)+1;
    % MEE0 = MEE; MEE(index,:) = a*MEE(index,:); 
    % param.MEE = MEE;
    % param.MEE_unperturbed = MEE0;
    % param.perturbation_type = 'local-rowwise(postsyn)';
    % param.perturbation_strength = a;    
    % MEE(:,index) = a*MEE(:,index); 
    % param.MEE = MEE;
    % param.perturbation_type = 'row+col perturbation';
    

%     MEE0 = MEE; MEE(:,index) = a*MEE(:,index); 
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'local-colwise(presyn)';
%     param.perturbation_strength = a;    
%     % % smooth local perturbation
%     index_x = 0.125*pi;
%     width_x = pi/4
%     perturbation = 1 - (1-a)*exp(-((x-index_x)/width_x).^2);
%     param.perturbation = perturbation;
%     perturbation = repmat(perturbation',1,nx);
%     MEE0 = MEE; MEE = MEE.*perturbation;
%     param.MEE = MEE;
%     param.MEE_unperturbed = MEE0;
%     param.perturbation_type = 'local-rowwise(postsyn)';
%     param.perturbation_strength = a;    
%     param.perturbation_index = index_x;
%     param.perturbation_width = width_x;
% global perturbation
    MEE0 = MEE; MEE = MEE*a;
    param.MEE = MEE;
    param.MEE_unperturbed = MEE0;
    param.perturbation_type = 'Global';
    param.perturbation_strength = a;
% random perturbation
%     a = 0.01;
%     param.perturbation_strength = a;
%     perturbation = 10.^(randn(nx)*a);
%     MEE0 = MEE; MEE = MEE.*perturbation;
%     param.perturbation_type = 'MEE-random-lognormal';
%     a = 0.03;
%     param.perturbation_strength = a;
%     perturbation = randn(nx)*a+1;
%     MEE0 = MEE; MEE = MEE.*perturbation;
%     param.perturbation_type = 'MEE-random-normal';
%     r = [0.9 1.1]; 
%     param.perturbation_range = r;
%     perturbation = rand(nx)*(r(2)-r(1))+r(1);
%     MEE0 = MEE; MEE = MEE.*perturbation;
%     param.perturbation_type = 'MEE-random-uniform';
    % a = 0.03;
    % param.perturbation_strength = a;
    % perturbation = gamrnd(1/a^2,a^2,nx,nx);
    % MEE0 = MEE; MEE = MEE.*perturbation;
    % param.perturbation_type = 'MEE-random-gamma';
%
    % param.perturbation = perturbation;
    % param.MEE = MEE;
    % param.MEE_unperturbed = MEE0;

    % previousResult = load("C:\Users\golde\Documents\Research\data\FR_Curr_ring_RK4_distractor_with_Plasticity\190806_11_11_LinearLargePerturb\results.mat");
    % MEE0 = MEE;MEE = previousResult.MEEt(:,:,end);
    % param.MEE = MEE;
    % param.MEE_unperturbed = MEE0;
    % % All Perturb
%     a = 0.02;
%     param.perturbation_strength = a;
%     perturbation = gamrnd(1/a^2,a^2,nx,nx);param.perturbationMEE = perturbation; MEE0 = MEE; MEE = MEE.*perturbation;
%     perturbation = gamrnd(1/a^2,a^2,nx,nx);param.perturbationMEI = perturbation; MEI0 = MEI; MEI = MEI.*perturbation;
%     perturbation = gamrnd(1/a^2,a^2,nx,nx);param.perturbationMIE = perturbation; MIE0 = MIE; MIE = MIE.*perturbation;
%     perturbation = gamrnd(1/a^2,a^2,nx,nx);param.perturbationMII = perturbation; MII0 = MII; MII = MII.*perturbation;
%     param.perturbation_type = 'all-random-gamma';
%     %
%     
%     param.MEE = MEE;param.MEE_unperturbed = MEE0;
%     param.MEI = MEI;param.MEI_unperturbed = MEI0;
%     param.MIE = MIE;param.MIE_unperturbed = MIE0;
%     param.MII = MII;param.MII_unperturbed = MII0;
    %% External Input
    JEO = 2*J;
    IEO_init = 1.35*(exp(-(x/(pi/4)).^2)'+1*ones(nx,1));
    
    param.IEo = JEO*IEO_init;
    param.IIo = 0;
   
    %% simulation timing in milisecond

    T_on = 500; %ms
    Tstim = 500;
    Tmemory = 3000;
    Tforget = 0;
    
    param.nTrialBatch = 200;
    nTrialMax=nTrialMax; % number of trails
    tTrial = T_on+Tstim+Tmemory+Tforget; % length of a trial
    tMax = nTrialMax*tTrial;

    TStimOn = T_on:tTrial:tMax;

    param.nTrial = nTrialMax;
    param.TStimOn   = TStimOn;
    param.TrialOn   = TStimOn-T_on;
    param.TStimOff  = TStimOn+Tstim;
    param.TDelayOff = TStimOn+Tstim+Tmemory;
    param.TForgetOff= TStimOn+Tstim+Tmemory+Tforget;
    param.Tmax = tMax;
    param.tTrial = tTrial;
    param.dt_store = 100;

    %% randomize stimlus location
    stimLoc = randi(floor(nx/np),np,nTrialMax); % random location in each group (1-4)
    stimLoc = stimLoc + (0:floor(nx/np):(nx-1))'; % add level to each group
    stimLoc = stimLoc - nx/2; % center to zero
    stimLoc_theta = stimLoc/nx*2*pi;

    pNp = randi(np,nTrialMax);

    param.stimLoc = stimLoc;
    param.stimLoc_theta = stimLoc_theta;
    param.pNp = pNp;

    %% additional parameters for plasticity
    % x: nx by 1, x: post-syn, x': pre-syn
    param.fM_expr = '@(x,dx) ( -lrD .* dx ) * x'' '; %differential plasticity within trial
    param.fM = eval(param.fM_expr);
    param.LearningRateDifferential = lrD;
    
    param.LearningRateHomeostatic = lrH;
    param.r_target = r_target;
    
    
