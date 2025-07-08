%Preset input for prepare-and-measure BB84
% The overlap between key and test basis can be tuned in parameters
% Channel model including loss, and unitary rotation (bitflip, phase flip
% angles)

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = OurProt_asymptotic()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handles of description/channel/error-correction files (functions)
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('OurProt_LossyDescription');
    % channelModel=str2func('BB84LossyNoisyChannel');
    % channelModel=str2func('OurProt_LinDispChannel');
    % channelModel=str2func('OurProt_QuadDispChannel');
    %channelModel=str2func('OurProt_MixedQuadDispChannel');
    %channelModel=str2func('OurProt_FBSChannel');
    % channelModel=str2func('OurProt_MixedFBSChannel');
    %channelModel=str2func('OurProt_LossMixedFBSChannel');
    %channelModel=str2func('OurProt_LossyChannel');
    %channelModel=str2func('OurProt_NonCollectiveLinDispChannel');
    channelModel=str2func('OurProt_NonCollectiveQuadDispChannel');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    %parameters.names = ["loss", "num_pts", "f","fullstat", "a", "theta", "phi"];
    parameters.names = ["loss", "num_pts", "f","fullstat", "a", "disp_const"];

    % %Encoding parameters:
    % Omega0 = -0.0095;
    % Omega1 = 0.0095;
    % sigma_w = 0.0011;
    % tau0 = -110;
    % tau1 = 110;
    % sigma_t = 30;


    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    %A = (0:0.04:2);
    %parameters.scan.loss = 10*(0:1:100)/100; % total loss
    %parameters.scan.loss = 10*(0:1:80)/100; % total loss
    %parameters.scan.theta = (0:0.04:2)*pi;
    %parameters.scan.phi = (0:0.04:2)*pi;
    %parameters.scan.theta = [A(50)-0.00001, A(50), A(50)+0.00001]*pi;
    %parameters.scan.phi = [A(47)-0.00001, A(47), A(47)+0.00001]*pi;
    %parameters.scan.theta = (A(138)+0.00001)*pi;
    %parameters.scan.theta = (0:0.1:1)*pi;
    % parameters.scan.phi = D(21); %[A(155)-0.000001,A(155),A(155)+0.000001]*pi;
    %parameters.scan.num_pts = (1:5:51);
    %parameters.scan.eta0 = B(113); %(0:0.008:1); % mode-dependent transmissivity
    %parameters.scan.delta = C(14);%(0:0.004:0.5)*pi;
    %parameters.scan.disp_const = (0:0.5:100);
    parameters.scan.disp_const = (0:5:2000);
    % parameters.scan.disp_const = (1100:5:1300);

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    parameters.fixed.loss = 0; % total loss
    parameters.fixed.fullstat = 1; %using full statistics or using QBER/Gain observables only
    parameters.fixed.f = 1; %useless error correction efficiency parameter
    parameters.fixed.a = sqrt(1/2); % overlap between key and test
    %parameters.fixed.disp_const = 0;
    %parameters.fixed.ftover = ftoverlap(Omega0, sigma_w, tau0, sigma_t);
    parameters.fixed.num_pts=1;
    %parameters.fixed.phi = 0;
    %parameters.fixed.theta=0;

    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    
%     parameters.optimize.pz = [0.1,0.5,0.9];
%     parameters.optimize.f = [1.0,1.2,2.0];
end

%set the running options for the solvers
function solverOptions=setOptions()

    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration
    %2.output at each solver 1 FW iteration
    %3.show all SDP solver output
    solverOptions.globalSetting.verboseLevel = 2; 
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
%     solverOptions.optimizer.name = 'localSearch_Adam';
    solverOptions.optimizer.linearResolution = 3; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 1; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 2; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 5; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 2; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.name = 'asymptotic';
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 1e-6; %1e-6 for asymptotic, 2.5e-3 for finite;
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'rref'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'asymptotic';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end