%Preset input for prepare-and-measure BB84
% The overlap between key and test basis can be tuned in parameters
% Channel model including loss, and unitary rotation (bitflip, phase flip
% angles)

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = BB84_asymptotic()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handles of description/channel/error-correction files (functions)
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('BB84LossyDescription');
    % channelModel=str2func('BB84LossyNoisyChannel');
    channelModel=str2func('BB84Channel1');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    parameters.names = ["loss", "eta0", "delta", "f","fullstat", "a", "theta", "phi"];
    %parameters.names = ["loss", "eta0", "f","fullstat", "a", "num_pts", "phi"];
    % parameters.names = ["loss", "num_pts", "f","fullstat", "a"];

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    
    % %A = (0:0.01:2);
    % B = (0:0.008:1);
    % C = (0:0.004:0.5)*pi;
    % D = (0:0.004:0.5)*pi;

    parameters.scan.loss = 10*(0:10:700)/100; % total loss
    % parameters.scan.theta = D(108); %[A(139)-0.000001,A(139),A(139)+0.000001]*pi;
    % parameters.scan.phi = D(21); %[A(155)-0.000001,A(155),A(155)+0.000001]*pi;
    % parameters.scan.num_pts = (0:100:1000);
    %parameters.scan.eta0 = B(113); %(0:0.008:1); % mode-dependent transmissivity
    %parameters.scan.delta = C(14);%(0:0.004:0.5)*pi;

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    %parameters.fixed.loss = 0; % total loss
    parameters.fixed.fullstat = 1; %using full statistics or using QBER/Gain observables only
    parameters.fixed.f = 1; %useless error correction efficiency parameter
    parameters.fixed.eta0 = 1.0; % mode-dependent transmissivity
    parameters.fixed.delta = 0.0;
    parameters.fixed.a = sqrt(1/2); % overlap between key and test
    parameters.fixed.theta = 0.0; % noise angle causing bitflip
    parameters.fixed.phi = 0.0; % noise angle causing phase flip

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