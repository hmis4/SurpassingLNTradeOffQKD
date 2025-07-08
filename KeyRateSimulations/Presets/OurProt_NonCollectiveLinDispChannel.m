%% FUNCTION NAME: BB84LossyNoisyChannel
% Lossy and noisy channel model for prepare-and-measure BB84. 
% Parameters:
% Transmitance eta, bitflip angle theta, phaseflip angle phi
%%

function channelModel = OurProt_NonCollectiveLinDispChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    %varNames=["loss","fullstat", "a", "num_pts", "ftover", "theta", "phi"];
    varNames = ["loss", "num_pts", "fullstat", "a", "disp_const"];
    
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
    dimPB = 5;%6; % Z basis (2), X basis (2), rest (1), no click (1)
    eta = 10^(-0.1*loss);  % transmissivity
    b = sqrt(1-a^2);

    %Encoding parameters:
    Omega0 = 0;
    Omega1 = 2*pi*0.019;
    sigma_w = 2*pi*0.0011;
    tau0 = 0;
    tau1 = 220;
    sigma_t = 30;

    delta_alpha = 5;

    
    [alpha_vec, beta_vec] = getOrigOverlaps(Omega0, Omega1,sigma_w,tau0,tau1,sigma_t);
    [ket0f, ket1f, ket0t, ket1t] = getKets(alpha_vec, beta_vec);

    alpha_00 = alpha_vec(1);
    alpha_01 = alpha_vec(2);
    alpha_10 = alpha_vec(3);
    alpha_11 = alpha_vec(4);

    ol_2d = alpha_00*alpha_11 - alpha_10*alpha_01;

    %Logical states
    Z0 = (kron(ket0f, ket1f)-kron(ket1f, ket0f))/sqrt(2);
    Z1 = ((kron(ket0t, ket1t)-kron(ket1t, ket0t))/sqrt(2) - ol_2d*Z0)/(sqrt(1-abs(ol_2d)^2));


    %Alice prepares
    signalStates = {Z0, Z1, (a*Z0 + b*Z1), (b*Z0-a*Z1)};  % 4 states
    probList = [1/4; 1/4; 1/4; 1/4];
    
    % rho_A constraints
    rhoA = zeros(dimA);
    %partial trace over flying qubit system to obtain local rhoA
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList(jRow) * probList(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
        end
    end
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(basis{iBasisElm}' * rhoA),'mask',0);
    end
    
    % Normalization
    addExpectations(1,'mask',0);
    
    rho_in = 0.25*kron(diag([1, 0, 0, 0]), blkdiag(signalStates{1}*signalStates{1}', [0])) + ...
    0.25*kron(diag([0, 1, 0, 0]),blkdiag(signalStates{2}*signalStates{2}', [0])) + ...
    0.25*kron(diag([0, 0, 1, 0]),blkdiag(signalStates{3}*signalStates{3}', [0])) + ...
    0.25*kron(diag([0, 0, 0, 1]),blkdiag(signalStates{4}*signalStates{4}', [0]));


    %%% LINEAR DISPERSION
    kappa_00 = probff1(disp_const, Omega0, Omega0, sigma_w);
    kappa_01 = probff1(disp_const, Omega0, Omega1, sigma_w);
    kappa_10 = probff1(disp_const, Omega1, Omega0, sigma_w);
    kappa_11 = probff1(disp_const, Omega1, Omega1, sigma_w);

    lambda_00 = probtt1(disp_const, tau0, tau0, sigma_t);
    lambda_01 = probtt1(disp_const, tau0, tau1, sigma_t);
    lambda_10 = probtt1(disp_const, tau1, tau0, sigma_t);
    lambda_11 = probtt1(disp_const, tau1, tau1, sigma_t);

    gamma_00 = probtf1(disp_const, Omega0, sigma_w, -tau0, sigma_t); %use minus sign to convert delta_ij to gamma_ij
    gamma_01 = probtf1(disp_const, Omega0, sigma_w, -tau1, sigma_t);
    gamma_10 = probtf1(disp_const, Omega1, sigma_w, -tau0, sigma_t);
    gamma_11 = probtf1(disp_const, Omega1, sigma_w, -tau1, sigma_t);
    
    delta_00 = probtf1(disp_const, Omega0, sigma_w, tau0, sigma_t); %designed as delta_ij (plus sign)
    delta_01 = probtf1(disp_const, Omega1, sigma_w, tau0, sigma_t);
    delta_10 = probtf1(disp_const, Omega0, sigma_w, tau1, sigma_t);
    delta_11 = probtf1(disp_const, Omega1, sigma_w, tau1, sigma_t);

    gamma_vec = [gamma_00, gamma_01, gamma_10, gamma_11];
    delta_vec = [delta_00, delta_01, delta_10, delta_11];
    kappa_vec = [kappa_00, kappa_01, kappa_10, kappa_11];
    lambda_vec = [lambda_00, lambda_01, lambda_10, lambda_11];

    F1 = getFop(gamma_vec, delta_vec, kappa_vec, lambda_vec, alpha_vec, beta_vec);


    kappa_00 = probff1(disp_const+delta_alpha, Omega0, Omega0, sigma_w);
    kappa_01 = probff1(disp_const+delta_alpha, Omega0, Omega1, sigma_w);
    kappa_10 = probff1(disp_const+delta_alpha, Omega1, Omega0, sigma_w);
    kappa_11 = probff1(disp_const+delta_alpha, Omega1, Omega1, sigma_w);

    lambda_00 = probtt1(disp_const+delta_alpha, tau0, tau0, sigma_t);
    lambda_01 = probtt1(disp_const+delta_alpha, tau0, tau1, sigma_t);
    lambda_10 = probtt1(disp_const+delta_alpha, tau1, tau0, sigma_t);
    lambda_11 = probtt1(disp_const+delta_alpha, tau1, tau1, sigma_t);

    gamma_00 = probtf1(disp_const+delta_alpha, Omega0, sigma_w, -tau0, sigma_t); %use minus sign to convert delta_ij to gamma_ij
    gamma_01 = probtf1(disp_const+delta_alpha, Omega0, sigma_w, -tau1, sigma_t);
    gamma_10 = probtf1(disp_const+delta_alpha, Omega1, sigma_w, -tau0, sigma_t);
    gamma_11 = probtf1(disp_const+delta_alpha, Omega1, sigma_w, -tau1, sigma_t);
    
    delta_00 = probtf1(disp_const+delta_alpha, Omega0, sigma_w, tau0, sigma_t); %designed as delta_ij (plus sign)
    delta_01 = probtf1(disp_const+delta_alpha, Omega1, sigma_w, tau0, sigma_t);
    delta_10 = probtf1(disp_const+delta_alpha, Omega0, sigma_w, tau1, sigma_t);
    delta_11 = probtf1(disp_const+delta_alpha, Omega1, sigma_w, tau1, sigma_t);

    gamma_vec = [gamma_00, gamma_01, gamma_10, gamma_11];
    delta_vec = [delta_00, delta_01, delta_10, delta_11];
    kappa_vec = [kappa_00, kappa_01, kappa_10, kappa_11];
    lambda_vec = [lambda_00, lambda_01, lambda_10, lambda_11];

    F2 = getFop(gamma_vec, delta_vec, kappa_vec, lambda_vec, alpha_vec, beta_vec);


    rho0 = kron(F1,F2)*(signalStates{1}*signalStates{1}')*(kron(F1,F2))';
    rho1 = kron(F1,F2)*(signalStates{2}*signalStates{2}')*(kron(F1,F2))';
    rhoplus = kron(F1,F2)*(signalStates{3}*signalStates{3}')*(kron(F1,F2))';
    rhominus = kron(F1,F2)*(signalStates{4}*signalStates{4}')*(kron(F1,F2))';

    rho_out = 0.25*kron(diag([1, 0, 0, 0]), blkdiag(rho0, [1-trace(rho0)])) + ...
    0.25*kron(diag([0, 1, 0, 0]),blkdiag(rho1, [1-trace(rho1)])) + ...
    0.25*kron(diag([0, 0, 1, 0]),blkdiag(rhoplus, [1-trace(rhoplus)])) + ...
    0.25*kron(diag([0, 0, 0, 1]),blkdiag(rhominus, [1-trace(rhominus)]));


    %Bob's POVMs
    basicBobPOVMs = {0.5*blkdiag(signalStates{1}*signalStates{1}', [0]), ...
        0.5*blkdiag(signalStates{2}*signalStates{2}', [0]),...
        0.5*blkdiag(signalStates{3}*signalStates{3}', [0]), ...
        0.5*blkdiag(signalStates{4}*signalStates{4}', [0]), ...
        blkdiag(eye(16) - 0.5*(signalStates{1}*signalStates{1}' + signalStates{2}*signalStates{2}' + signalStates{3}*signalStates{3}' + signalStates{4}*signalStates{4}'), [1]), ...
        %blkdiag(zeros(16), [1])
        };

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
    
    % 4 states
    gain=sum(Y1_simulation.*[1,1,0,0,0; 1,1,0,0,0; 0,0,1,1,0; 0,0,1,1,0],'all');
    error=sum(Y1_simulation.*[0,1,0,0,0; 1,0,0,0,0; 0,0,1,0,0; 0,0,0,1,0],'all')/gain;


    
    % my_proba_dist = Y1_simulation(1:2,1:2,1);
    % psift=sum(sum(my_proba_dist));
    % 
    % %channelModel.errorRate = [error]; % to compute EC based on qber
    % channelModel.probDist = my_proba_dist/psift; % to compute EC on dist


    my_proba_dist_z = Y1_simulation(1:2,1:2,1);
    psift_z=sum(sum(my_proba_dist_z));

    my_proba_dist_x = Y1_simulation(3:4,3:4,1);
    psift_x=sum(sum(my_proba_dist_x));

    %channelModel.errorRate = [error]; % to compute EC based on qber
    channelModel.probDist = {my_proba_dist_z/psift_z, my_proba_dist_x/psift_x}; % to compute EC on dist
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    
    % channelModel.pSift = [gain];  % to compute EC based on qber
    channelModel.pSift = [psift_z, psift_x];  % EC based on proba dist
end
