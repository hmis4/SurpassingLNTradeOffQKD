%% FUNCTION NAME: BoileauLossyChannel
% Lossy channel model for prepare-and-measure Boileau protocol 3 photons decoherence free
%%

function channelModel = Boileau3LossyChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["loss","fullstat"];
    
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
    dimPB = 10;

    eta = 10^(-0.1*loss);  % erasure
    
    % basic qubit povm Z
    Z0 = [1,0;0,0];
    Z1 = [0,0;0,1];

    % basic qubit povm X
    X0 = [1,1;1,1]/2;
    X1 = [1,-1;-1,1]/2;

    % list all measurement operators according to Hannah's notes
    % they do not correspond to projectors on z or x states but rather pinched projector
    Pi_az = kron(kron(Z0, Z0), Z0) + kron(kron(Z1, Z1), Z1);
    Pi_ax = kron(kron(X0, X0), X0) + kron(kron(X1, X1), X1);

    Pi_bz = kron(kron(Z0, Z0), Z1) + kron(kron(Z1, Z1), Z0);
    Pi_bx = kron(kron(X0, X0), X1) + kron(kron(X1, X1), X0);

    Pi_cz = kron(kron(Z0, Z1), Z0) + kron(kron(Z1, Z0), Z1);
    Pi_cx = kron(kron(X0, X1), X0) + kron(kron(X1, X0), X1);

    Pi_dz = kron(kron(Z0, Z1), Z1) + kron(kron(Z1, Z0), Z0);
    Pi_dx = kron(kron(X0, X1), X1) + kron(kron(X1, X0), X0);

    % list all 3 qubit states
    z000 = kron(kron(zket(2,1), zket(2,1)), zket(2,1));
    z001 = kron(kron(zket(2,1), zket(2,1)), zket(2,2));
    z010 = kron(kron(zket(2,1), zket(2,2)), zket(2,1));
    z011 = kron(kron(zket(2,1), zket(2,2)), zket(2,2));
    z100 = kron(kron(zket(2,2), zket(2,1)), zket(2,1));
    z101 = kron(kron(zket(2,2), zket(2,1)), zket(2,2));
    z110 = kron(kron(zket(2,2), zket(2,2)), zket(2,1));
    z111 = kron(kron(zket(2,2), zket(2,2)), zket(2,2));


    signalStates = {
        (z100-z010)/sqrt(2),...
        (z101-z011)/sqrt(2),...
        (z100-z001)/sqrt(2), ...
        (z110-z011)/sqrt(2), ...
        (z100-z001)/sqrt(2), ...
        (z110-z011)/sqrt(2), ...
        (z010-z001)/sqrt(2), ...
        (z110-z101)/sqrt(2)};  % 4 states, middle one is repeated, all are duplicated for purification
    probList = ones(1,8)/8;
    
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
    rho_in = 0.125*kron(diag([1,0,0,0,0,0,0,0]), blkdiag(signalStates{1}*signalStates{1}', [0])) + ...
        0.125*kron(diag([0,1,0,0,0,0,0,0]),blkdiag(signalStates{2}*signalStates{2}', [0])) + ...
        0.125*kron(diag([0,0,1,0,0,0,0,0]),blkdiag(signalStates{3}*signalStates{3}', [0])) + ...
        0.125*kron(diag([0,0,0,1,0,0,0,0]),blkdiag(signalStates{4}*signalStates{4}', [0])) + ...
        0.125*kron(diag([0,0,0,0,1,0,0,0]),blkdiag(signalStates{5}*signalStates{5}', [0])) + ...
        0.125*kron(diag([0,0,0,0,0,1,0,0]),blkdiag(signalStates{6}*signalStates{6}', [0])) + ...
        0.125*kron(diag([0,0,0,0,0,0,1,0]),blkdiag(signalStates{7}*signalStates{7}', [0])) + ...
        0.125*kron(diag([0,0,0,0,0,0,0,1]),blkdiag(signalStates{8}*signalStates{8}', [0]));

    rho_out = eta^3*rho_in + (1-eta^3)*kron(eye(8)/8, blkdiag(zeros(8), [1])); % 3 photon protocol
    

    basicBobPOVMs = {0.5*(blkdiag(Pi_az, [0])), ...
        0.5*(blkdiag(Pi_bz, [0])),...
        0.5*(blkdiag(Pi_cz, [0])), ...
        0.5*(blkdiag(Pi_dz, [0])), ...
        0.5*(blkdiag(Pi_ax, [0])), ...
        0.5*(blkdiag(Pi_bx, [0])), ...
        0.5*(blkdiag(Pi_cx, [0])), ...
        0.5*(blkdiag(Pi_dx, [0])), ...
        blkdiag(eye(8) - 0.5*(Pi_az + Pi_bz + Pi_cz + Pi_dz + Pi_ax + Pi_bx + Pi_cx + Pi_dx), [0]), ...
        blkdiag(zeros(8), [1]), 
        };  % all POVMs

    basicAlicePOVMs = {
        diag([1,0,0,0,0,0,0,0]),...
        diag([0,1,0,0,0,0,0,0]), ...
        diag([0,0,1,0,0,0,0,0]), ...
        diag([0,0,0,1,0,0,0,0]),...
        diag([0,0,0,0,1,0,0,0]),...
        diag([0,0,0,0,0,1,0,0]),...
        diag([0,0,0,0,0,0,1,0]),...
        diag([0,0,0,0,0,0,0,1])
        }; % 4 states, each duplicated

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
    
    % 4 states + duplicates
    % for key basis, mean az or ax for bit 0 and cz or cx for bit 1
    % old version
    % gain=sum(Y1_simulation.*[1,0,1,0,1,0,1,0,0; 1,0,1,0,1,0,1,0,0; ...
    %     1,0,1,0,1,0,1,0,0; 1,0,1,0,1,0,1,0,0; ...
    %     0,0,0,0,0,0,0,0,0; 0,0,0,0,0,0,0,0,0; ...
    %     0,0,0,0,0,0,0,0,0; 0,0,0,0,0,0,0,0,0],'all');
    % error=sum(Y1_simulation.*[0,0,1,0,0,0,1,0,0;0,0,1,0,0,0,1,0,0; ...
    %     1,0,0,0,1,0,0,0,0;1,0,0,0,1,0,0,0,0; ...
    %     0,0,0,0,0,0,0,0,0; 0,0,0,0,0,0,0,0,0; ...
    %     0,0,0,0,0,0,0,0,0; 0,0,0,0,0,0,0,0,0],'all')/gain;

    stats = reshape(sum(reshape(Y1_simulation, 2,4, 10), 1), 4, 10);

    my_proba_dist_z = [stats(1:2,2:3,1),stats(1:2,6:7,1)];
    psift_z=sum(sum(my_proba_dist_z));

    my_proba_dist_x = [stats(3:4,3:4,1),stats(3:4,7:8,1)];
    psift_x=sum(sum(my_proba_dist_x));


    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    % channelModel.errorRate = [error];
    % channelModel.pSift = [gain];
    channelModel.pSift = [psift_z, psift_x];
    channelModel.probDist = {my_proba_dist_z/psift_z, my_proba_dist_x/psift_x}; % to compute EC on dist

end
