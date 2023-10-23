function pZ = quadMap(pZ,varargin)
% quadMap - computes the quadratic map of a polyZonotope
%
% Syntax:
%    pZ = quadMap(pZ,Q)
%    pZ = quadMap(pZ,pZ2,Q)
%    pZ = quadMap(pZ,pZ2,Q,dep)
%
% Inputs:
%    pZ, pZ2 - polyZonotope objects
%    Q - quadratic coefficients as a cell of matrices
%    dep - keep dependencies (dep = true) or not (dep = false) 
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    % quadratic multiplication
%    pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);
%    
%    Q{1} = [1 2; -1 2];
%    Q{2} = [-3 0; 1 1];
%
%    pZquad = quadMap(pZ,Q);
%
%    figure;
%    subplot(1,2,1);
%    plot(pZ,[1,2],'FaceColor','r');
%    subplot(1,2,2);
%    plot(pZquad,[1,2],'FaceColor','b');
%
%    % mixed quadratic multiplication
%    pZ2 = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[0 0 0;0 0 0;1 0 3;0 1 1]);
%
%    pZquadMixed = quadMap(pZ,pZ2,Q);
%
%    figure
%    subplot(1,3,1);
%    plot(pZ,[1,2],'FaceColor','r');
%    subplot(1,3,2);
%    plot(pZ2,[1,2],'FaceColor','b');
%    subplot(1,3,3);
%    plot(pZquadMixed,[1,2],'FaceColor','g','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/quadMap, cubMap

% Authors:       Niklas Kochdumper
% Written:       23-March-2018
% Last update:   21-April-2020 (remove zero-length independent generators)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if nargin == 1
        throw(CORAerror('CORAerror:notEnoughInputArgs',2));
    elseif nargin == 2
        pZ = aux_quadMapSingle(pZ,varargin{1});
    elseif nargin == 3
        pZ = aux_quadMapMixed(pZ,varargin{1},varargin{2},false);
    elseif nargin == 4
        if ~isscalar(varargin{3}) || ~islogical(varargin{3})
            throw(CORAerror('CORA:wrongValue', 'fourth', ...
                                            'has to be boolean.'));
        end
        pZ = aux_quadMapMixed(pZ,varargin{1},varargin{2},varargin{3});
    else
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end
end


% Auxiliary functions -----------------------------------------------------

function pZ = aux_quadMapSingle(pZ,Q)
% compute an over-approximation of the quadratic map 
%
% {x_i = x^T Q{i} x | x \in pZ} 

    % split into a zonotope Z that overapproximates the dependent generators,
    % and a zonotope Zrem that contains the independent generators
    % (Z + Zrem)' Q (Z + Zrem) = ...
    %       Z' Q Z  +  Zrem' Q Z  +  Z' Q Zrem  +  Zrem' Q Zrem 
    
    pZtemp = pZ;
    pZtemp.GI = [];

    Z = zonotope(pZtemp);
    Zrem = zonotope([0*pZ.c, pZ.GI]);
    
    % return directly if pZ is a zonotope (rest matrix only)
    if isempty(pZ.G)
        Zrem_ = zonotope([pZ.c, pZ.GI]);
        Zr = quadMap(Zrem_,Q);
        pZ.c = Zr.c;
        pZ.GI = Zr.G;
        return;
    end

    % extend generator and exponent matrix by center
    Gext = [pZ.c, pZ.G];
    Eext = [zeros(size(pZ.E,1),1), pZ.E];


    % initialize the resulting generator and exponent matrix 
    N = size(Gext,2);
    dim = length(Q);
    M = N*(N+1)/2;

    Equad = zeros(size(pZ.E,1),M);
    Gquad = zeros(dim,M);

    % create the exponent matrix that corresponds to the quadratic map
    counter = 1;

    for j = 1:N
        Equad(:,counter:counter+N-j) = Eext(:,j:N) + Eext(:,j)*ones(1,N-j+1);
        counter = counter + N - j + 1;
    end


    % loop over all dimensions
    for i = 1:dim
        if ~isempty(Q{i})
            % quadratic evaluation
            quadMat = Gext'*Q{i}*Gext;
    
            % copy the result into the generator matrix
            counter = 1;
    
            for j = 1:N
                Gquad(i,counter:counter+N-j) = ...
                    [quadMat(j,j), quadMat(j,j+1:N) + quadMat(j+1:N,j)'];
                counter = counter + N - j + 1;
            end
        end
    end

    % add up all generators that belong to identical exponents
    [Enew,Gnew] = removeRedundantExponents(Equad,Gquad);
    
    
    % quadratic and mixed multiplication of remaining generators
    if ~isempty(pZ.GI)

        % quadratic multiplication
        Ztemp1 = quadMap(Zrem,Q);

        % mixed multiplications
        Ztemp2 = quadMap(Z,Zrem,Q);
        Ztemp3 = quadMap(Zrem,Z,Q);

        pZ.c = Ztemp1.c + Ztemp2.c + Ztemp3.c;
        GI = [Ztemp1.G, Ztemp2.G, Ztemp3.G];
        
        % delete generators of length zero
        GI = GI(:,any(GI,1));

    else
        pZ.c = zeros(length(Q),1);
        GI = [];
    end

    % assemble the properties of the resulting polynomial zonotope
    if sum(Enew(:,1)) == 0
        pZ.c = pZ.c + Gnew(:,1);
        pZ.G = Gnew(:,2:end);
        pZ.E = Enew(:,2:end);
        pZ.GI = GI;
    else
        pZ.G = Gnew;
        pZ.E = Enew;
        pZ.GI = GI;
    end
end

function pZ = aux_quadMapMixed(pZ1,pZ2,Q,dep)
% compute an over-approximation of the quadratic map 
%
% {x_i = x1^T Q{i} x2 | x1 \in pZ1, x2 \in pZ2} 
%
% of two polyZonotope objects.

    % bring the exponent matrices to a common representation
    if ~dep
        pZ2.id = max(pZ1.id) + pZ2.id;
    end

    [id,E1,E2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.E,pZ2.E);
    id = (1:length(id))';
    
    % split into a zonotope Z that overapproximates the dependent generators,
    % and a zonotope Zrem that contains the independent generators
    pZtemp = pZ1;
    pZtemp.GI = [];
    Z1 = zonotope(pZtemp);
    Zrem1 = zonotope([0*pZ1.c, pZ1.GI]);
    
    pZtemp = pZ2;
    pZtemp.GI = [];
    Z2 = zonotope(pZtemp);
    Zrem2 = zonotope([0*pZ2.c, pZ2.GI]);

    % construct extended generator and exponent matrix (extendet by center)
    Gext1 = [pZ1.c, pZ1.G];
    Eext1 = [zeros(size(E1,1),1), E1];
    
    Gext2 = [pZ2.c, pZ2.G];
    Eext2 = [zeros(size(E2,1),1), E2];

    % initialize the resulting generator and exponent matrix 
    N1 = size(Gext1,2);
    N2 = size(Gext2,2);
    
    dim = length(Q);
    M = N1*N2;

    Equad = zeros(size(E1,1),M);
    Gquad = zeros(dim,M);

    % create the exponent matrix that corresponds to the quadratic map
    counter = 1;

    for j = 1:N2
        Equad(:,counter:counter+N1-1) = Eext1 + Eext2(:,j)*ones(1,N1);
        counter = counter + N1 ;
    end
    
    % loop over all dimensions
    for i = 1:length(Q)
        if ~isempty(Q{i})
            % quadratic evaluation
            quadMat = Gext1'*Q{i}*Gext2;
            Gquad(i,:) = reshape(quadMat,1,[]);
        end
    end

    % add up all generators that belong to identical exponents
    [Enew,Gnew] = removeRedundantExponents(Equad,Gquad);

    % mixed multiplication of remaining generators
    GI = [];
    Zquad_rest = quadMap(Zrem1,Zrem2,Q);
    G_add = Zquad_rest.G;
    c = Zquad_rest.c;
    GI(:,end+1:end+size(G_add,2)) = G_add;

    Zquad_mixed1 = quadMap(Z1,Zrem2,Q);
    G_add = Zquad_mixed1.G;
    c = c + Zquad_rest.c;
    GI(:,end+1:end+size(G_add,2)) = G_add;

    Zquad_mixed2 = quadMap(Zrem1,Z2,Q);
    G_add = Zquad_mixed2.G;
    c = c + Zquad_rest.c;
    GI(:,end+1:end+size(G_add,2)) = G_add;
    
    % remove zero-length independent generators
    GI = GI(:,any(GI,1));

    % assemble the properties of the resulting polynomial zonotope
    if sum(Enew(:,1)) == 0
        pZ = polyZonotope(c + Gnew(:,1),Gnew(:,2:end),GI,Enew(:,2:end),id);
    else
        pZ = polyZonotope(c,Gnew,GI,Enew,id);
    end

end

% ------------------------------ END OF CODE ------------------------------
