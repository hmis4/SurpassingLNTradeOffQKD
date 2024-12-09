
%% FUNCTION NAME: Boileau3LossyDescription
% Lossy description for prepare-and-measure Boileau protocol, 3 photons, decoherence free
%%

function protocolDescription = Boileau3LossyDescription(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["fullstat"];
    
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
    
    dimA = 8; % prepared mixed states, purification double the dimension
    dimB = 9;  % dim 8+1 (3 qubits+loss)
    dimPB = 10; % a,b,c,d in each basis and rest, noclick
    

    % define all states and measurements
    az0 = kron(kron(zket(2,1), zket(2,1)), kron(zket(2,1), zket(2,1)));
    az1 = kron(kron(zket(2,2), zket(2,2)), kron(zket(2,2), zket(2,2)));

    bz0 = kron(kron(zket(2,1), zket(2,1)), kron(zket(2,2), zket(2,2)));
    bz1 = kron(kron(zket(2,2), zket(2,2)), kron(zket(2,1), zket(2,1)));

    cz0 = kron(kron(zket(2,1), zket(2,2)), kron(zket(2,1), zket(2,2)));
    cz1 = kron(kron(zket(2,2), zket(2,1)), kron(zket(2,2), zket(2,1)));

    dz0 = kron(kron(zket(2,1), zket(2,2)), kron(zket(2,2), zket(2,1)));
    dz1 = kron(kron(zket(2,2), zket(2,1)), kron(zket(2,1), zket(2,2)));
    
    zStates = {(az0 + az1)/sqrt(2), ...
        (bz0+ bz1)/sqrt(2), ...
        (cz0 + cz1)/sqrt(2), ...
        (dz0 + dz1)/sqrt(2)
        }; % a,b,c,d in z basis
    
    H_gate = [1,1;1,-1]/sqrt(2);
    H4_gate = kron(kron(H_gate, H_gate),kron(H_gate, H_gate));
    ax0 = H4_gate * az0;
    ax1 = H4_gate * az1;

    bx0 = H4_gate * bz0;
    bx1 = H4_gate * bz1;

    cx0 = H4_gate * cz0;
    cx1 = H4_gate * cz1;

    dx0 = H4_gate * dz0;
    dx1 = H4_gate * dz1;

    xStates = {(ax0 + ax1)/sqrt(2), ...
        (bx0+ bx1)/sqrt(2), ...
        (cx0 + cx1)/sqrt(2), ...
        (dx0 + dx1)/sqrt(2)
        }; % a,b,c,d in x basis
    
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

    
    
    % 2 key states are duplicated because of purification
    % old version key from one basis only
    % krausOpC = kron(kron(kron(zket(2,1),diag([1,1,0,0,0,0,0,0])) ...
    %     + kron(zket(2,2),diag([0,0,1,1,0,0,0,0])), ...
    %     blkdiag(sqrtBobconclusivePovm,[0])), [1]);  % conclusive
    % krausOp = {krausOpC};
 
    % Boileau protocol, keymap on Bob. 
    % z basis alice, zbasis bob
    krausOpCZZ = kron(kron(kron(zket(2,2),diag([1,1,1,1,0,0,0,0])), blkdiag(Pi_bz/sqrt(2),[0])), zket(4,1))+...
        kron(kron(kron(zket(2,1),diag([1,1,1,1,0,0,0,0])), blkdiag(Pi_cz/sqrt(2),[0])), zket(4,1));
    
    % x basis alice, zbasis bob
    krausOpCXZ = kron(kron(kron(zket(2,2),diag([0,0,0,0,1,1,1,1])), blkdiag(Pi_cz/sqrt(2),[0])), zket(4,3))+...
        kron(kron(kron(zket(2,1),diag([0,0,0,0,1,1,1,1])), blkdiag(Pi_dz/sqrt(2),[0])), zket(4,3));
    
    % z basis alice, xbasis bob
    krausOpCZX = kron(kron(kron(zket(2,2),diag([1,1,1,1,0,0,0,0])), blkdiag(Pi_bx/sqrt(2),[0])), zket(4,2))+...
        kron(kron(kron(zket(2,1),diag([1,1,1,1,0,0,0,0])), blkdiag(Pi_cx/sqrt(2),[0])), zket(4,2));

    % x basis alice, xbasis bob
    krausOpCXX = kron(kron(kron(zket(2,2),diag([0,0,0,0,1,1,1,1])), blkdiag(Pi_cx/sqrt(2),[0])), zket(4,4))+...
        kron(kron(kron(zket(2,1),diag([0,0,0,0,1,1,1,1])), blkdiag(Pi_dx/sqrt(2),[0])), zket(4,4));
    krausOp = {krausOpCZZ, krausOpCZX, krausOpCXZ, krausOpCXX};
 
    % components for the pinching Z map
    % one basis
    % keyProj1 = kron(diag([1, 0]), eye(dimA*dimB)); 
    % keyProj2 = kron(diag([0, 1]), eye(dimA*dimB));
    
    % 2 basis
    keyProj1 = kron(diag([1, 0]), eye(4*dimA*dimB)); 
    keyProj2 = kron(diag([0, 1]), eye(4*dimA*dimB));

    keyMap = {keyProj1, keyProj2};

    % Constraints
    %observables = {};

    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addObservables(kron(basis{iBasisElm}, eye(dimB)),'mask',0);
    end
    
    % Normalization
    addObservables(eye(dimA*dimB),'mask',0);
    
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