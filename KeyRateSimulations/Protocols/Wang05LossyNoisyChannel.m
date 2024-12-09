%% FUNCTION NAME: Wang05LossyNoisyChannel
% Lossy channel model for prepare-and-measure Wang05 protocol. 
%%

function channelModel = Wang05LossyNoisyChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["loss","fullstat", "theta", "phi"];
    %varNames=["loss","fullstat", "num_pts", "phi"];
    
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
    dimPB = 6;
    eta = 10^(-0.1*loss);  % transmissivity
    signalStates = {[0;1;0;0], [0;0;1;0], [0;1;1;0]/sqrt(2), [0;1;-1;0]/sqrt(2)};  % 4 states
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

    % unitary noise
    U1 = [cos(theta), -sin(theta)*exp(-1i*phi); 
        sin(theta)*exp(1i*phi), cos(theta)];

    U = kron(diag([1,1,1,1]),blkdiag(kron(U1, U1), [1]));  % collective noise on Bob's system
    
    % lossy channel
    rho_in = 0.25*kron(diag([1, 0, 0, 0]), diag([0,1,0,0,0])) ...
    + 0.25*kron(diag([0,1,0,0]), diag([0,0,1,0,0])) ...
    + 0.25*kron(diag([0,0,1,0]), [0,0,0,0,0;0,1,1,0,0;0,1,1,0,0;0,0,0,0,0;0,0,0,0,0]/2) ...
    + 0.25*kron(diag([0,0,0,1]), [0,0,0,0,0;0,1,-1,0,0;0,-1,1,0,0;0,0,0,0,0;0,0,0,0,0]/2) ; % 4 states

    rho_out = eta^2*U*rho_in*(U') + (1-eta^2)*kron(diag([1,1,1,1])/4,diag([0,0,0,0,1])); % eta^2 because 2 photon protocol
    
    % % Mixed channel
    % theta_vec = linspace(0,2*pi,num_pts+1);
    % p_theta = 1/num_pts;
    % 
    % rho_out = (1-eta^2)*kron(diag([1,1,1,1])/4,diag([0,0,0,0,1]));
    % 
    % for i=1:num_pts
    %     U1 = [cos(theta_vec(i)), -sin(theta_vec(i))*exp(-1i*phi);
    %         sin(theta_vec(i))*exp(1i*phi), cos(theta_vec(i))];
    % 
    %     U = kron(diag([1,1,1,1]),blkdiag(kron(U1, U1), [1]));
    % 
    %     rho_out = rho_out + eta^2*(U*rho_in*U')*p_theta;
    % end

    basicBobPOVMs = {0.5*diag([0,1,0,0,0]), ...
        0.5*diag([0,0,1,0,0]),...
        0.5*0.5*[0,0,0,0,0;0,1,1,0,0;0,1,1,0,0;0,0,0,0,0;0,0,0,0,0], ...
        0.5*0.5*[0,0,0,0,0;0,1,-1,0,0;0,-1,1,0,0;0,0,0,0,0;0,0,0,0,0],...
        diag([1,0,0,1,0]), ...
        diag([0,0,0,0,1])};
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
    
    % qber but actually not use, rather use whole proba distrib below
    % gain=sum(Y1_simulation.*[1,1,0,0,0; 1,1,0,0,0; 0,0,1,1,0; 0,0,1,1,0],'all');
    % error=sum(Y1_simulation.*[0,1,0,0,0; 1,0,0,0,0; 0,0,0,1,0; 0,0,1,0,0],'all')/gain;


    my_proba_dist_z = Y1_simulation(1:2,1:2,1);
    psift_z=sum(sum(my_proba_dist_z));

    my_proba_dist_x = Y1_simulation(3:4,3:4,1);
    psift_x=sum(sum(my_proba_dist_x));
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    % channelModel.errorRate = [error]; % for qber computation
    % channelModel.pSift = [gain]; % for qber computation
    channelModel.pSift = [psift_z, psift_x];
    channelModel.probDist = {my_proba_dist_z/psift_z, my_proba_dist_x/psift_x}; % to compute EC on dist
    

end
