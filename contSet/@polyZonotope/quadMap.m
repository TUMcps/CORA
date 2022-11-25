function res = quadMap(varargin)
% quadMap - computes the quadratic map of a polyZonotope
%
% Syntax:  
%    res = quadMap(pZ,Q)
%    res = quadMap(pZ1,pZ2,Q)
%
% Inputs:
%    pZ,pZ1,pZ2 - polyZonotope objects
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    res - resulting set as a polyZonotope object
%
% Example: 
%    % quadrtic multiplication
%    pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);
%    
%    Q{1} = [1 2; -1 2];
%    Q{2} = [-3 0; 1 1];
%
%    pZquad = quadMap(pZ,Q);
%
%    figure
%    subplot(1,2,1);
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    subplot(1,2,2);
%    plot(pZquad,[1,2],'b','Filled',true,'EdgeColor','none');
%
%    % mixed quadratic multiplication
%    pZ2 = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[0 0 0;0 0 0;1 0 3;0 1 1]);
%
%    pZquadMixed = quadMap(pZ,pZ2,Q);
%
%    figure
%    subplot(1,3,1);
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    subplot(1,3,2);
%    plot(pZ2,[1,2],'b','Filled',true,'EdgeColor','none');
%    subplot(1,3,3);
%    plot(pZquadMixed,[1,2],'g','Filled',true,'EdgeColor','none','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/quadMap, cubMap

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  21-April-2020 (remove zero-length independent generators)
% Last revision:---

%------------- BEGIN CODE --------------

    if nargin == 2
        res = quadMapSingle(varargin{1},varargin{2});
    else
        res = quadMapMixed(varargin{1},varargin{2},varargin{3});
    end
end


% Auxiliary Functions -----------------------------------------------------

function pZ = quadMapSingle(pZ,Q)
% compute an over-approximation of the quadratic map 
%
% {x_i = x^T Q{i} x | x \in pZ} 

    % split into a zonotope Z that overapproximates the dependent generators,
    % and a zonotope Zrem that contains the independent generators
    % (Z + Zrem)' Q (Z + Zrem) = ...
    %       Z' Q Z  +  Zrem' Q Z  +  Z' Q Zrem  +  Zrem' Q Zrem 
    pZtemp = pZ;
    pZtemp.Grest = [];

    Z = zonotope(pZtemp);
    Zrem = zonotope([0*pZ.c, pZ.Grest]);

    % extend generator and exponent matrix by center
    Gext = [pZ.c, pZ.G];
    Eext = [zeros(size(pZ.expMat,1),1), pZ.expMat];


    % initialize the resulting generator and exponent matrix 
    N = size(Gext,2);
    dim = length(Q);
    M = N*(N+1)/2;

    Equad = zeros(size(pZ.expMat,1),M);
    Gquad = zeros(dim,M);

    % create the exponent matrix that corresponds to the quadratic map
    counter = 1;

    for j = 1:N
        Equad(:,counter:counter+N-j) = Eext(:,j:N) + Eext(:,j)*ones(1,N-j+1);
        counter = counter + N - j + 1;
    end


    % loop over all dimensions
    for i = 1:dim

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

    % add up all generators that belong to identical exponents
    [ExpNew,Gnew] = removeRedundantExponents(Equad,Gquad);
    
    
    % quadratic and mixed multiplication of remaining generators
    if ~isempty(pZ.Grest)

        % quadratic multiplication
        Ztemp1 = quadMap(Zrem,Q);

        % mixed multiplications
        Ztemp2 = quadMap(Z,Zrem,Q);
        Ztemp3 = quadMap(Zrem,Z,Q);

        pZ.c = center(Ztemp1) + center(Ztemp2) + center(Ztemp3);
        Grest = [generators(Ztemp1), generators(Ztemp2), generators(Ztemp3)];
        
        % delete generators of length zero
        Grest = Grest(:,any(Grest,1));

    else
        pZ.c = zeros(length(Q),1);
        Grest = [];
    end

    % assemble the properties of the resulting polynomial zonotope
    if sum(ExpNew(:,1)) == 0
        pZ.c = pZ.c + Gnew(:,1);
        pZ.G = Gnew(:,2:end);
        pZ.expMat = ExpNew(:,2:end);
        pZ.Grest = Grest;
    else
        pZ.G = Gnew;
        pZ.expMat = ExpNew;
        pZ.Grest = Grest;
    end
end

function pZ = quadMapMixed(pZ1,pZ2,Q)
% compute an over-approximation of the quadratic map 
%
% {x_i = x1^T Q{i} x2 | x1 \in pZ1, x2 \in pZ2} 
%
% of two polyZonotope objects.

    % bring the exponent matrices to a common representation
    pZ2.id = max(pZ1.id) + pZ2.id;
    [id,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);
    id = (1:length(id))';
    
    % split into a zonotope Z that overapproximates the dependent generators,
    % and a zonotope Zrem that contains the independent generators
    pZtemp = pZ1;
    pZtemp.Grest = [];
    Z1 = zonotope(pZtemp);
    Zrem1 = zonotope([0*pZ1.c, pZ1.Grest]);
    
    pZtemp = pZ2;
    pZtemp.Grest = [];
    Z2 = zonotope(pZtemp);
    Zrem2 = zonotope([0*pZ2.c, pZ2.Grest]);

    % construct extended generator and exponent matrix (extendet by center)
    Gext1 = [pZ1.c, pZ1.G];
    Eext1 = [zeros(size(expMat1,1),1), expMat1];
    
    Gext2 = [pZ2.c, pZ2.G];
    Eext2 = [zeros(size(expMat2,1),1), expMat2];

    % initialize the resulting generator and exponent matrix 
    N1 = size(Gext1,2);
    N2 = size(Gext2,2);
    
    dim = length(Q);
    M = N1*N2;

    Equad = zeros(size(expMat1,1),M);
    Gquad = zeros(dim,M);

    % create the exponent matrix that corresponds to the quadratic map
    counter = 1;

    for j = 1:N2
        Equad(:,counter:counter+N1-1) = Eext1 + Eext2(:,j)*ones(1,N1);
        counter = counter + N1 ;
    end
    
    % loop over all dimensions
    for i = 1:length(Q)
           
        % quadratic evaluation
        quadMat = Gext1'*Q{i}*Gext2;
        Gquad(i,:) = reshape(quadMat,1,[]);        
    end

    % add up all generators that belong to identical exponents
    [ExpNew,Gnew] = removeRedundantExponents(Equad,Gquad);

    % mixed multiplication of remaining generators
    Grest = [];
    Zquad_rest = quadMap(Zrem1,Zrem2,Q);
    G_add = generators(Zquad_rest);
    c = center(Zquad_rest);
    Grest(:,end+1:end+length(G_add(1,:))) = G_add;

    Zquad_mixed1 = quadMap(Z1,Zrem2,Q);
    G_add = generators(Zquad_mixed1);
    c = c + center(Zquad_rest);
    Grest(:,end+1:end+length(G_add(1,:))) = G_add;

    Zquad_mixed2 = quadMap(Zrem1,Z2,Q);
    G_add = generators(Zquad_mixed2);
    c = c + center(Zquad_rest);
    Grest(:,end+1:end+length(G_add(1,:))) = G_add;
    
    % remove zero-length generators
    Grest = Grest(:,any(Grest,1));

    % assemble the properties of the resulting polynomial zonotope
    if sum(ExpNew(:,1)) == 0
        pZ = polyZonotope(c + Gnew(:,1),Gnew(:,2:end),Grest,ExpNew(:,2:end),id);
    else
        pZ = polyZonotope(c,Gnew,Grest,ExpNew,id);
    end

end

%------------- END OF CODE --------------