
%% FUNCTION NAME: pmBB84LossyDescription
% Lossy description for prepare-and-measure BB84.
% The observables correspond to a squashing model with five POVM outcomes
% (including photon loss).
% Detectors are assumed perfect (no dark count, no misalignment, no loss)
% Any actual loss or misalignment is attributed to the adversary
%%

function protocolDescription = OurProt_LossyDescription(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["fullstat", "a"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    observables = {};
    obsMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % dimA = 2; % 2 states
    dimA = 4; % 4 states
    dimB = 17;  % 2 ququarts + loss
    dimPB = 5;%6; % Z basis (2), X basis (2), rest (1), no click (1)
    b = sqrt(1-a^2);


    %Encoding parameters:
    Omega0 = 0;
    Omega1 = 2*pi*0.019;
    sigma_w = 2*pi*0.0011;
    tau0 = 0;
    tau1 = 220;
    sigma_t = 30;
    
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

    X0 = (a*Z0 + b*Z1);
    X1 = (b*Z0-a*Z1);
    %signalStates = {Z0, Z1, (a*Z0 + b*Z1), (b*Z0-a*Z1)};  % 4 states

    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
    % they include postselection into the Kraus operators directly


    krausOpCZ = kron(kron(kron(zket(2,1),diag([1,1,0,0])), blkdiag(Z0*Z0',[0])/sqrt(2)) + ...
                    kron(kron(zket(2,2),diag([1,1,0,0])), blkdiag(Z1*Z1', [0])/sqrt(2)), [1;0]);  % conclusive Z

    krausOpCX = kron(kron(kron(zket(2,1),diag([0,0,1,1])), blkdiag(X0*X0', [0])/sqrt(2)) + ...
                    kron(kron(zket(2,2),diag([0,0,1,1])), blkdiag(X1*X1', [0])/sqrt(2)), [0;1]);  % conclusive X

    krausOp = {krausOpCZ, krausOpCX};
 
    % key from 2 bases
    keyProj1 = kron(diag([1, 0]), eye(dimA*dimB*2)); 
    keyProj2 = kron(diag([0, 1]), eye(dimA*dimB*2));
    keyMap = {keyProj1, keyProj2};

    % Constraints
    %observables = {};

    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addObservables(kron(basis{iBasisElm}, eye(dimB)),'mask',0);
    end
    
    % Normalization
    addObservables(eye(dimA*dimB),'mask',0);
    
    %Constructing the constraints
    % basicBobPOVMs = {0.5*diag([1,0,0]), ...
    %     0.5*diag([0,1,0]),...
    %     0.5*[a^2,a*b,0;a*b,b^2,0;0,0,0], ...
    %     0.5*[b^2,-a*b,0;-a*b,a^2,0;0,0,0], ...
    %     diag([0,0,1])};

    basicBobPOVMs = {0.5*blkdiag(Z0*Z0', [0]), ...
        0.5*blkdiag(Z1*Z1', [0]),...
        0.5*blkdiag(X0*X0', [0]), ...
        0.5*blkdiag(X1*X1', [0]), ...
        blkdiag(eye(16) - 0.5*(Z0*Z0' + Z1*Z1' + X0*X0' + X1*X1'), [1]), ...
        %blkdiag(zeros(16), [1]),
        };

    basicAlicePOVMs = {diag([1,0,0,0]),diag([0,1,0,0]), diag([0,0,1,0]), diag([0,0,0,1])}; % 4 states
    
    %Full set of bipartite POVMS
    bipartitePOVMs = cell(dimA*dimPB,1);
    for i = 1:dimA
        for j = 1:dimPB
            bipartitePOVMs{dimPB*(i-1)+(j-1)+1} = kron(basicAlicePOVMs{i},basicBobPOVMs{j});
        end
    end    
    
    if(fullstat==1)
        addObservables(bipartitePOVMs,'mask',1);
    else
        %QBER and Gain statistics
        select=@(x,y)dimPB*(x-1)+(y-1)+1;
        %observables = [observables; bipartitePOVMs{select(1,2)} ; bipartitePOVMs{select(2,1)}; bipartitePOVMs{select(3,4)} ; bipartitePOVMs{select(4,3)}];
        newObservables = {bipartitePOVMs{select(1,2)} ; bipartitePOVMs{select(2,1)}; bipartitePOVMs{select(1,1)} ; bipartitePOVMs{select(2,2)}; ...
            bipartitePOVMs{select(3,4)} ; bipartitePOVMs{select(4,3)}; bipartitePOVMs{select(3,3)} ; bipartitePOVMs{select(4,4)}};
        addObservables(newObservables,'mask',1);
        
        %normalization
        temp = bipartitePOVMs{select(1,2)}+bipartitePOVMs{select(2,1)}+bipartitePOVMs{select(1,1)}+bipartitePOVMs{select(2,2)}+ ...
            bipartitePOVMs{select(3,4)}+bipartitePOVMs{select(4,3)}+ bipartitePOVMs{select(3,3)} + bipartitePOVMs{select(4,4)};
        
        addObservables(eye(dimA*dimB)-temp,'mask',1);
    end
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];

end