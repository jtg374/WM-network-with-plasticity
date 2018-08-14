
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
    nx = 80;
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

    %% Perturbations
    a = .85;
    MEE0 = MEE; MEE = a*MEE; 
    param.MEE = MEE;
    param.MEE_unperturbed = MEE0;
    param.perturbation_type = 'global';
    param.perturbation_strength = a;
    
    
    %% External Input
    JEO = 2*J;
    IEO_init = 1.35*(exp(-(x/(pi/4)).^2)'+1*ones(nx,1));
    
    param.IEo = JEO*IEO_init;
    param.IIo = 0;
    param.JWipe=2.5;
    
    %% simulation timing in milisecond

    T_on = 500; %ms
    Tstim = 500;
    Tmemory = 3000;
    Tforget = Tstim*2;
    iter=500;
    
    Tmax = T_on+iter*(Tstim+Tmemory+Tforget)-Tforget; % end of training
    Tinit = T_on:(Tstim+Tmemory+Tforget):Tmax; % Times of stimuli onset (training peroid) 

    Tinter = 6000;  % memory time in test period
    Tmax1 = Tmax+Tforget; % Stimilus onset time in test period    
    Tmax2 = Tmax1+Tstim+Tinter; 
    
    param.nTrailTrain = iter;
    param.nTrailTest = 1;
    param.TStimOn   = [Tinit, Tmax1];
    param.TStimOff  = [Tinit+Tstim, Tmax1+Tstim];
    param.TDelayOff = [Tinit+Tstim+Tmemory, Tmax2];
    param.TForgetOff= param.TStimOn(2:end);
    param.Tmax = Tmax2;

    param.dt_store = 50;

    %% randomize stimlus location
    stimLoc = randi(nx,iter,1)-nx/2; % in training period
    stimLoc_theta = stimLoc/nx*2*pi;
    stimLoc_test = 0; 

    param.stimLoc = [stimLoc; stimLoc_test];
    param.stimLoc_theta = [stimLoc_theta; stimLoc_test/nx*2*pi];

    %% additional parameters for plasticity
    % x: nx by 1, x: post-syn, x': pre-syn
    param.fM_expr = "@(x,dx) ( ((x/20).^(2*0.5)+10) .* dx ) * x'";
    param.fM = eval(param.fM_expr);
