function param = NDF_with_Plasticity_Parameters(a,b,lrD,lrH,nTrial,r_target)
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
    JEE = 1; JEI = 1; JIE = 1;
    sigma_E = 0.2*pi; sigma_I = 1*pi; % wider excitatory synaptic projection
    fEE = exp(-(x/sigma_E).^2)/(2*pi);
    fIE = exp(-(x/sigma_E).^2)/(2*pi);
    fEI = exp(-(x/sigma_I).^2)/(2*pi);


    MEE = zeros(nx,nx);
    MIE = zeros(nx,nx);
    MEI = zeros(nx,nx);

    for i = 1:nx
        MEE(i,:) = JEE*dx*circshift(fEE,[0 -round(pi/(dx))-1+i]);
        MIE(i,:) = JIE*dx*circshift(fIE,[0 -round(pi/(dx))-1+i]);
        MEI(i,:) = JEI*dx*circshift(fEI,[0 -round(pi/(dx))-1+i]);
    end

    param.MEE = MEE;
    param.MIE = MIE;
    param.MEI = MEI;

    %% Transfer Function
    qE = @(x) 7*((0.2*(x-1)+0.5).*(x<=1) + (2*(x-1)+0.5).*(x>1).*(x<=2) + (1*(x-2)+2.5).*(x>2).*(x<14) + 14.5.*(x>14));
    qI = @(x) x;
    param.qE = qE;
    param.qI = qI;

    %% Perturbations
%     % % smooth local perturbation
    center_x = 0*pi;
    width_x = pi/4;
    perturbation_row = 1 - (1-a)*exp(-((x-center_x)/width_x).^2);
    perturbation_col = 1 - (1-b)*exp(-((x-center_x)/width_x).^2);
    perturbation = perturbation_row'*perturbation_col;
    param.perturbation = perturbation;
    MEE0 = MEE; MEE = MEE.*perturbation;
    param.MEE = MEE;
    param.MEE_unperturbed = MEE0;
    param.perturbation_type = 'local-smooth';
    param.perturbation_strength_row = a;    
    param.perturbation_strength_col = b;    
    param.perturbation_index = center_x;
    param.perturbation_width = width_x;
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
    param.IEc = 0.6*ones(nx,1);
    
    JEO = 1;
    sigma_o = 0.2*pi;
    IEO_init = 0.5*(exp(-(x/sigma_o).^2)'+1*ones(nx,1));
    
    param.IEo = JEO*IEO_init;
    param.IIo = 0;
   
    %% simulation timing in milisecond

    T_on = 500; %ms
    Tstim = 500;
    Tmemory = 3000;
    Tforget = 0;
    
    nTrialMax=nTrial; % number of trails
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
    param.fM_XS_expr = '@(x,dx)  -lrD .* dx*x'' '; %differential plasticity 
    param.fM_homeo_expr = '@(x,g)  lrH * (r_target-x).*g '; %homeostatic plasticity
    param.fM_XS = eval(param.fM_XS_expr);
    param.fM_homeo = eval(param.fM_homeo_expr);
    param.LearningRateDifferential = lrD;
    
    param.LearningRateHomeostatic = lrH;
    param.r_target = r_target;
    
    
