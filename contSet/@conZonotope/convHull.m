function cZ = convHull(cZ,S,varargin)
% convHull - computes the convex hull of a constrained zonotope and one or
%    multiple other set representations
%
% Syntax:
%    cZ = convHull(cZ,S)
%    cZ = convHull(cZ,{S1,...,Sm})
%    cZ = convHull(cZ,S,method)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object (or cell-array of conZonotope objects)
%    method - (optional) method for computation of convex hull
%                        'exact:Raghuraman' (default)
%                        'exact:Kochdumper'
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
% 
%    Z = [4 2 0 0;4 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ2 = conZonotope(Z,A,b);
%
%    res = convHull(cZ1,cZ2);
%
%    figure; hold on;
%    plot(cZ1,[1,2],'FaceColor','r');
%    plot(cZ2,[1,2],'FaceColor','b');
%    plot(res,[1,2],'g','LineWidth',3);
%
% Reference:
%    [1] V. Raghuraman, J.P. Koeln, 'Set operations and order reductions
%        for constrained zonotopes', arXiv:2009.06039v1
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/enclose

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       13-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
%                21-December-2022 (MW, implemented method [Thm.5,1], change syntax for conZonotope objects)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % only one set given
    if nargin == 1
        return
    end

    % too many input arguments
    if nargin > 3
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end

    % set default values
    method = setDefaultValues({'exact:Raghuraman'},varargin);

    % use different algorithms for the case with only one or two sets
    if ~iscell(S)

        if representsa_(S,'emptySet',eps)
            return
        elseif representsa_(cZ,'emptySet',1e-10)
            cZ = S; return
        end

        % find a conZonotope object
        [cZ,S] = findClassArg(cZ,S,'conZonotope');

        % handle different classes of the second set
        if isa(S,'conZonotope')

            if strcmp(method,'exact:Raghuraman')
                cZ = aux_convHullcZ_Raghuraman(cZ,S);
            elseif strcmp(method,'exact:Kochdumper')
                cZ = aux_convHullcZ_Kochdumper(cZ,S);
            end

        elseif isa(S,'zonotope') || isa(S,'interval') || ...
               isa(S,'polytope') || isa(S,'zonoBundle') || isnumeric(S)

            if strcmp(method,'exact:Raghuraman')
                cZ = aux_convHullcZ_Raghuraman(cZ,conZonotope(S));
            elseif strcmp(method,'exact:Kochdumper')
                cZ = aux_convHullcZ_Kochdumper(cZ,conZonotope(S));
            end
            
        elseif isa(S,'polyZonotope') || isa(S,'conPolyZono')
            
            % convert conZonotope to polyZonotope
            cZ = convHull(polyZonotope(cZ),S);

        else
            % throw error for given arguments
            throw(CORAerror('CORA:noops',cZ,S));

        end

    else % S is a cell-array
        
        if all(cellfun(@(x) isa(x,'conZonotope'),S,'UniformOutput',true))

            % multiple sets only if all are constrained zonotopes
            cZ = aux_convHullMany(cZ,S);

        else
            % throw error
            throw(CORAerror('CORA:noops',cZ,varargin(:)));

        end
    end
end


% Auxiliary functions -----------------------------------------------------

function cZ = aux_convHullcZ_Raghuraman(cZ1,cZ2)
% implement method from [Thm.5,1]

    % obtain centers and generator matrices of constrained zonotopes
    c1 = cZ1.c; G1 = cZ1.G;
    c2 = cZ2.c; G2 = cZ2.G;
    
    % obtain dimension and number of constraints
    [n,nrGens1] = size(G1);
    nrGens2 = size(G2,2);
    
    % obtain constraint matrices, offsets, and number of constraints
    A1 = cZ1.A; b1 = cZ1.b; nrCon1 = size(A1,1);
    A2 = cZ2.A; b2 = cZ2.b; nrCon2 = size(A2,1);

    % compute center of convex hull
    c = 0.5*(c1+c2);

    % compute generater matrix
    G = [G1, G2, 0.5*(c1-c2), zeros(n,2*(nrGens1+nrGens2))];

    % compute sparsity of resulting constraint matrix
    minNumZeros = nrCon1*nrGens2 + nrCon1*2*(nrGens1+nrGens2) ... first row in A
        + nrCon2*nrGens1 + nrCon2*2*(nrGens1+nrGens2)         ... second row in A
        + 2*nrGens1*(nrGens1-1) + 4*nrGens2*nrGens1 + 2*nrGens2*(nrGens2-1) ...
        + (2*(nrGens1+nrGens2))^2 - 2*(nrGens1+nrGens2);      ... third row in A
    numElem = (nrCon1 + nrCon2 + 2*nrGens1+2*nrGens2) ...
        * (nrGens1 + nrGens2 + 1 + 2*(nrGens1+nrGens2));
    sparsity = minNumZeros / numElem;

    % compute constraint matrix
    if n < 10 || sparsity < 0.5
        A31 = [eye(nrGens1); -eye(nrGens1); zeros(nrGens2,nrGens1); zeros(nrGens2,nrGens1)];
        A32 = [zeros(nrGens1,nrGens2); zeros(nrGens1,nrGens2); eye(nrGens2); -eye(nrGens2)];
        A30 = [-0.5*ones(nrGens1,1); -0.5*ones(nrGens1,1); 0.5*ones(nrGens2,1); 0.5*ones(nrGens2,1)];
        A = [A1, zeros(nrCon1,nrGens2), -0.5*b1, zeros(nrCon1,2*(nrGens1+nrGens2));
             zeros(nrCon2,nrGens1), A2, 0.5*b2, zeros(nrCon2,2*(nrGens1+nrGens2));
             A31, A32, A30, eye(2*(nrGens1+nrGens2))];

    else
        % use sparse matrices for large number of constraints/generators
        A31_sparse = [speye(nrGens1); -speye(nrGens1); sparse(nrGens2,nrGens1); sparse(nrGens2,nrGens1)];
        A32_sparse = [sparse(nrGens1,nrGens2); sparse(nrGens1,nrGens2); speye(nrGens2); -speye(nrGens2)];
        A30_sparse = [-0.5*ones(2*nrGens1,1); 0.5*ones(2*nrGens2,1)];
        A_part1 = [A1, sparse(nrCon1,nrGens2), -0.5*b1, sparse(nrCon1,2*(nrGens1+nrGens2))];
        A_part2 = [sparse(nrCon2,nrGens1), A2, 0.5*b2, sparse(nrCon2,2*(nrGens1+nrGens2));];
        A_part3 = [A31_sparse, A32_sparse, A30_sparse, speye(2*(nrGens1+nrGens2))];
        A = [A_part1; A_part2; A_part3];
    end

    % compute offset vector
    b = [0.5*b1; 0.5*b2; -0.5*ones(2*(nrGens1+nrGens2),1)];

    % instantiate resulting constrained zonotope
    cZ = conZonotope(c,G,A,b);

end

function cZ = aux_convHullcZ_Kochdumper(cZ1,cZ2)
% compute convex hull of two constrained zonotopes

    % obtain object properties
    c1 = cZ1.c; G1 = cZ1.G;
    c2 = cZ2.c; G2 = cZ2.G;

    % obtain dimension and number of constraints
    [n,m1] = size(G1);
    m2 = size(G2,2);
    
    % obtain constraint matrices (if given)
    if isempty(cZ1.A)
        A1 = zeros(1,m1); b1 = 0;
    else
        A1 = cZ1.A; b1 = cZ1.b;
    end
    
    if isempty(cZ2.A)
        A2 = zeros(1,m2); b2 = 0;
    else
        A2 = cZ2.A; b2 = cZ2.b;
    end
    
    % obtain number of constraints
    p1 = size(A1,1);
    p2 = size(A2,1);

    % compute center
    c = 0.5 * (c1 + c2);

    % compute generator matrix
    G = [G1 G2 0.5*(c1-c2) zeros(n,3*m1 + 3*m2)];

    % compute constraint matrix
    A1_ = [blkdiag(A1,A2),0.5*[-b1;b2]];
    A2_ = [-2*eye(2*m1), [-2*eye(m1);-2*eye(m1)]];
    A3_ = [-2*eye(2*m2), [-2*eye(m2);-2*eye(m2)]];
    
    A = blkdiag(A1_,A2_,A3_);
    
    A(p1+p2+1:p1+p2+2*m1,1:m1) = [-2*eye(m1);2*eye(m1)];
    A(p1+p2+2*m1+1:end,m1+1:m1+m2) = [2*eye(m2);-2*eye(m2)];
    A(p1+p2+1:end,m1+m2+1) = [ones(2*m1,1);-ones(2*m2,1)];

    % compute constraint offset
    b = [0.5*b1;0.5*b2;3*ones(2*m1+2*m2,1)];
    
    % remove trivial constraints 0 = 0
    ind = find(sum(abs([A,b]),2) > 0);
    A = A(ind,:);
    b = b(ind);

    % construct resulting constraint zonotope object
    cZ = conZonotope([c,G],A,b);
end


function cZ = aux_convHullMany(cZ1,list)
% compute convex hull of many constrained zonotopes

    % initialize variables
    n = dim(cZ1);
    
    A_ = cell(length(list),1);
    b_ = [];
    a = [];
    G_ = [];
    c_ = zeros(n,1);
    
    % loop over all sets
    for i = 1:length(list)
        
        % obtain object properties
        A = list{i}.A;
        b = list{i}.b;
        G = list{i}.G;
        c = list{i}.c;
        
        p = size(A,1);
        m = size(G,2);
        
        % construct constraint matrix and vector
        I = eye(m);
        O = zeros(m);
        O_ = zeros(p,m);
        o = ones(m,1);
        
        if ~isempty(A)
            A_{i} = [-0.5*b A O_ O_ O_;o -2*I -2*I O -2*I;o 2*I O -2*I -2*I];
            b_ = [b_;0.5*b;3*ones(2*m,1)];
        else
            A_{i} = [o -2*I -2*I O -2*I;o 2*I O -2*I -2*I];
            b_ = [b_;3*ones(2*m,1)];
        end
        
        a = [a, 0.5, zeros(1,4*m)];
        
        % construct generator matrix
        G_ = [G_,0.5*c,G,zeros(n,3*m)];
        
        % construct center 
        c_ = c_ + 0.5*c;
        
    end
    
    % construct the resulting conZonotope object
    A = [blkdiag(A_{:});a];
    b = [b_;1-0.5*length(list)];
    
    cZ = conZonotope([c_,G_],A,b);  

end

% ------------------------------ END OF CODE ------------------------------
