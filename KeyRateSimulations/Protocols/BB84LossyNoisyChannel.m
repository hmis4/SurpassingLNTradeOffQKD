%% FUNCTION NAME: BB84LossyNoisyChannel
% Lossy and noisy channel model for prepare-and-measure BB84. 
% Parameters:
% Transmitance eta, bitflip angle theta, phaseflip angle phi
%%

function channelModel = BB84LossyNoisyChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["loss","fullstat", "a", "theta", "phi"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this channel model file
    %can be automatically filled in by calling addExpectations(x) or addExpectations(x,'mask',maskValue)
    expectations = [];
    expMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model begin %%%%%%%%%%%%%%%%%%%%%%%%%

    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    dimPB = 5;
    eta = 10^(-0.1*loss);  % transmissivity
    b = sqrt(1-a^2);
    signalStates = {[1;0], [0; 1], [a;b], [b;-a]};  % 4 states
    probList = [1/4; 1/4; 1/4; 1/4];
    
    % rho_A constraints
    rhoA = zeros(dimA);
    %partial trace over flying qubit system to obtain local rhoA
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList(jRow) * probList(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
        end
    end
    expectations = [];
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(basis{iBasisElm}' * rhoA),'mask',0);
    end
    
    % Normalization
    addExpectations(1,'mask',0);
    
    % lossy channel
    rho_in = 0.25*kron(diag([1, 0, 0, 0]),diag([1, 0, 0])) ...
    + 0.25*kron(diag([0,1,0,0]), diag([0,1,0])) ...
    + 0.25*kron(diag([0,0,1,0]), [a^2,a*b,0 ;a*b,b^2,0; 0,0,0]) ...
    + 0.25*kron(diag([0,0,0,1]), [b^2,-a*b,0; -a*b,a^2,0; 0,0,0]) ; % 4 states
    
    % Noise unitary, misalignment
    U1 = [cos(theta), -sin(theta)*exp(1i*phi); 
        sin(theta), cos(theta)*exp(1i*phi)];

    U = kron(diag([1,1,1,1]),blkdiag(U1, [1]));

    % one photon protocol
    rho_out = eta^1*(U*rho_in*U') + (1-eta^1)*kron(diag([1,1,1,1])/4,diag([0,0,1]));
    
    basicBobPOVMs = {0.5*diag([1,0,0]), ...
        0.5*diag([0,1,0]),...
        0.5*[a^2,a*b,0;a*b,b^2,0;0,0,0], ...
        0.5*[b^2,-a*b,0;-a*b,a^2,0;0,0,0], ...
        diag([0,0,1])};
    basicAlicePOVMs = {diag([1,0,0,0]),diag([0,1,0,0]), diag([0,0,1,0]), diag([0,0,0,1])};
    %Full set of bipartite POVMS
    for i = 1:dimA
        for j = 1:dimPB
            POVM=kron(basicAlicePOVMs{i},basicBobPOVMs{j});
            Y1_simulation(i,j)=real(trace(POVM'*rho_out));
        end
    end
    

    bipartiteExpectations = Y1_simulation;  % stats from simulation
    bipartiteExpectations_1D = zeros(dimA*dimPB,1);
    for i = 1:dimA
        for j = 1:dimPB
            bipartiteExpectations_1D(dimPB*(i-1)+(j-1)+1) = bipartiteExpectations(i,j);
        end
    end
    
    
    if(fullstat==1)
        %fine-grain statistics
        addExpectations(bipartiteExpectations_1D,'mask',1);
    else
        %QBER and Gain statistics
        select=@(x,y)dimPB*(x-1)+(y-1)+1;
        %     expectations = [expectations; bipartiteExpectations_1D(select(1,2));bipartiteExpectations_1D(select(2,1));bipartiteExpectations_1D(select(3,4));bipartiteExpectations_1D(select(4,3))];
        coarseGrainExpectations = [bipartiteExpectations_1D(select(1,2));bipartiteExpectations_1D(select(2,1));bipartiteExpectations_1D(select(1,1));bipartiteExpectations_1D(select(2,2));...
        bipartiteExpectations_1D(select(3,4));bipartiteExpectations_1D(select(4,3));bipartiteExpectations_1D(select(3,3));bipartiteExpectations_1D(select(4,4))];
        addExpectations(coarseGrainExpectations,'mask',1);
        %normalization
        temp = sum([bipartiteExpectations_1D(select(1,2));bipartiteExpectations_1D(select(2,1));bipartiteExpectations_1D(select(1,1));bipartiteExpectations_1D(select(2,2));...
        bipartiteExpectations_1D(select(3,4));bipartiteExpectations_1D(select(4,3));bipartiteExpectations_1D(select(3,3));bipartiteExpectations_1D(select(4,4))]);
        addExpectations(1-temp,'mask',1);
    end
    
    % 2 states
    
    % 4 states
    gain=sum(Y1_simulation.*[1,1,0,0,0; 1,1,0,0,0; 0,0,0,0,0; 0,0,0,0,0],'all');
    error=sum(Y1_simulation.*[0,1,0,0,0; 1,0,0,0,0; 0,0,0,0,0; 0,0,0,0,0],'all')/gain;
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [error];
    channelModel.pSift = [gain];

end
