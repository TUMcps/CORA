function pZ = cubMap(pZ, varargin)
% cubMap - computes the set corresponding to the cubic multiplication of a
%    polynomial zonotope with a third-order tensor
%
% Description:
%    Calculates the following set:
%    { z = (x' T x) * x | x \in pZ }
%
%    If three polyZonotopes are provided, the function calculates the set:
%    { z = (x1' T x2) * x3 | x1 \in pZ1, x2 \in pZ2, x3 \in pZ3 }
%
% Syntax:
%    pZ = cubMap(pZ,T)
%    pZ = cubMap(pZ,T,ind)
%    pZ = cubMap(pZ1,pZ2,pZ3,T)
%    pZ = cubMap(pZ1,pZ2,pZ3,T,ind)
%
% Inputs:
%    pZ,pZ1,pZ2,pZ3 - polyZonotope objects
%    T - third-order tensor
%    ind - cell-array containing the non-zero indices of the tensor
%
% Outputs:
%    pZ - polyZonotope object representing the set of the cubic mapping
%
% Example:
%    % cubic multiplication
%    pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);
%
%    T{1,1} = [1 2; -1 2];
%    T{1,2} = [-3 0; 1 1];
%    T{2,1} = [2 0; -2 1];
%    T{2,2} = [-3 0; -21 -1];
%
%    pZcub = cubMap(pZ,T);
%
%    figure;
%    subplot(1,2,1)
%    plot(pZ,[1,2],'FaceColor','r');
%    subplot(1,2,2)
%    plot(pZcub,[1,2],'FaceColor','b');
%
%    % mixed cubic multiplication
%    pZ2 = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[0 0 0;0 0 0;1 0 3;0 1 1]);
%
%    pZcubMixed = cubMap(pZ,pZ2,pZ2,T);
%
%    figure
%    subplot(1,3,1);
%    plot(pZ,[1,2],'FaceColor','r');
%    subplot(1,3,2);
%    plot(pZ2,[1,2],'FaceColor','b');
%    subplot(1,3,3);
%    plot(pZcubMixed,[1,2],'FaceColor','g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadMap, zonotope/cubMap

% Authors:       Niklas Kochdumper
% Written:       17-August-2018
% Last update:   25-May-2023 (TL, bugfix unique ids for GI)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs', 2));
elseif nargin > 5
    throw(CORAerror('CORA:tooManyInputArgs', 5));
end

% cubic multiplication or mixed cubic multiplication
if nargin == 4 || nargin == 5
    % syntax cases:
    % res = cubMap(pZ1,pZ2,pZ3,T)
    % res = cubMap(pZ1,pZ2,pZ3,T,ind)

    % assign input arguments
    pZ2 = varargin{1};
    pZ3 = varargin{2};
    T = varargin{3};

    % parse optional input arguments
    if nargin == 5
        ind = varargin{4};
    else
        temp = 1:size(T, 2);
        ind = repmat({temp}, [size(T, 1), 1]);
    end

    % check input arguments
    inputArgsCheck({{pZ, 'att', 'polyZonotope'}; ...
        {pZ2, 'att', 'polyZonotope'}; ...
        {pZ3, 'att', 'polyZonotope'}; ...
        {T, 'att', 'cell'}; ...
        {ind, 'att', 'cell'}});

    % mixed cubic multiplication
    pZ = aux_cubMapMixed(pZ, pZ2, pZ3, T, ind);

elseif nargin == 2 || nargin == 3
    % syntax cases:
    % res = cubMap(pZ,T)
    % res = cubMap(pZ,T,ind)

    % assign input argument
    T = varargin{1};

    % parse optional input arguments
    if nargin == 3
        ind = varargin{2};
    else
        temp = 1:size(T, 2);
        ind = repmat({temp}, [size(T, 1), 1]);
    end

    % check input arguments
    inputArgsCheck({{pZ, 'att', 'polyZonotope'}; ...
        {T, 'att', 'cell'}; ...
        {ind, 'att', 'cell'}});

    % cubic multiplication
    pZ = aux_cubMapSingle(pZ, T, ind);
end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_cubMapSingle(pZ, T, ind)
% calculates the following set:      { z = (x' T x) * x | x \in pZ }

% split pZ into pZdep that contains the dependent generators,
% and pZind that contains the independent generators

% remove independent generators for pZdep
pZdep = pZ;
pZdep.GI = [];

% create new pZ with unique ids for indepedent generators
% this makes the computation more accurate but will be removed later
m = size(pZ.GI, 2);
maxId = max(pZ.id);
indIds = (maxId + 1:maxId + m)';
pZind = polyZonotope(0*pZ.c, pZ.GI, [], eye(m), indIds);

% construct extended generator and exponent matrix (extended by center)
Gext = [pZ.c, pZ.G];
Eext = [zeros(size(pZ.E, 1), 1), pZ.E];

% initialize the resulting generator and exponent matrix
N = size(Gext, 2);
n = length(ind);
M = N * (N + 1) / 2;

Equad = zeros(size(pZ.E, 1), M);
Ecub = zeros(size(Equad, 1), size(Equad, 2)*N);
Gcub = zeros(n, size(Equad, 2)*N);

% create the exponent matrix that corresponds to the quadratic map
counter = 1;
for j = 1:N
    Equad(:, counter:counter+N-j) = Eext(:, j:N) + Eext(:, j) * ones(1, N-j+1);
    counter = counter + N - j + 1;
end

% create the exponent matrix that corresponds to the cubic map
for j = 1:N
    Ecub(:, (j - 1)*M+1:j*M) = Equad + Eext(:, j) * ones(1, M);
end

% loop over all dimensions
for i = 1:length(ind)

    % initialize quadratic matrix
    Gquad = zeros(1, M);

    % loop over all quadratic matrices: \sum_k (pZ' T_k pZ) * pZ_k
    for k = 1:length(ind{i})

        % quadratic evaluation
        quadMat = Gext' * T{i, ind{i}(k)} * Gext;

        % copy the result into the generator matrix
        counter = 1;
        for j = 1:N
            Gquad(counter:counter+N-j) = ...
                [quadMat(j, j), quadMat(j, j+1:N) + quadMat(j+1:N, j)'];
            counter = counter + N - j + 1;
        end

        % cubic generator matrix: loop over all generators
        for j = 1:N
            Gcub(i, (j - 1)*M+1:j*M) = Gcub(i, (j - 1)*M+1:j*M) + ...
                Gquad * Gext(ind{i}(k), j);
        end
    end
end

% add up all generators that belong to identical exponents
[ExpNew, Gnew] = removeRedundantExponents(Ecub, Gcub);

% mixed multiplication with the zonotope from the independent terms
if ~isempty(pZ.GI) && ~all(all(pZ.GI == 0))

    % cubic map of pZind with itself
    pZindMap = cubMap(pZind, T, ind);

    % cubic map with all combinations between pZind and pZdep
    pZMapList = cell(6, 1);
    pZMapList{1} = cubMap(pZind, pZind, pZdep, T, ind);
    pZMapList{2} = cubMap(pZind, pZdep, pZdep, T, ind);
    pZMapList{3} = cubMap(pZind, pZdep, pZind, T, ind);
    pZMapList{4} = cubMap(pZdep, pZind, pZdep, T, ind);
    pZMapList{5} = cubMap(pZdep, pZind, pZind, T, ind);
    pZMapList{6} = cubMap(pZdep, pZdep, pZind, T, ind);

    % sum of all sets
    pZindSum = sum(pZindMap, pZMapList);

    % remove dependencies of previously independent generators
    Zind = zonotope(pZindSum);

    % get resulting c and GI
    GInew = Zind.G;
    cnew = Zind.c;

else
    GInew = [];
    cnew = zeros(n, 1);
end

% construct the resulting polynomial zonotope
if all(ExpNew(:, 1) == 0)
    res = polyZonotope(cnew+Gnew(:, 1), ...
        Gnew(:, 2:end), GInew, ExpNew(:, 2:end), pZ.id);
else
    res = polyZonotope(cnew, Gnew, GInew, ExpNew, pZ.id);
end
end

function res = aux_cubMapMixed(pZ1, pZ2, pZ3, T, ind)
% calculates the following set:
% { z = (x1' T x2) * x3 | x1 \in pZ1, x2 \in pZ2, x3 \in pZ3 }

% bring the exponent matrices to a common representation
[id_, E1, E2] = mergeExpMatrix(pZ1.id, pZ2.id, pZ1.E, pZ2.E);
[id, E1, E3] = mergeExpMatrix(id_, pZ3.id, E1, pZ3.E);
[~, E2, ~] = mergeExpMatrix(id_, pZ3.id, E2, pZ3.E);

% split into a zonotope Z that over-approximates the dependent generators,
% and a zonotope Zrem that contains the independent generators
pZtemp = pZ1;
pZtemp.GI = [];
Z1 = zonotope(pZtemp);
Zrem1 = zonotope([0 * pZ1.c, pZ1.GI]);

pZtemp = pZ2;
pZtemp.GI = [];
Z2 = zonotope(pZtemp);
Zrem2 = zonotope([0 * pZ2.c, pZ2.GI]);

pZtemp = pZ3;
pZtemp.GI = [];
Z3 = zonotope(pZtemp);
Zrem3 = zonotope([0 * pZ3.c, pZ3.GI]);

% construct extended generator and exponent matrix (extended by center)
Gext1 = [pZ1.c, pZ1.G];
Eext1 = [zeros(size(E1, 1), 1), E1];

Gext2 = [pZ2.c, pZ2.G];
Eext2 = [zeros(size(E2, 1), 1), E2];

Gext3 = [pZ3.c, pZ3.G];
Eext3 = [zeros(size(E3, 1), 1), E3];

% initialize the resulting generator and exponent matrix
N1 = size(Gext1, 2);
N2 = size(Gext2, 2);
N3 = size(Gext3, 2);

n = length(ind);
M = N1 * N2;

Equad = zeros(size(E1, 1), M);
Ecub = zeros(size(Equad, 1), size(Equad, 2)*N3);
Gcub = zeros(n, size(Equad, 2)*N3);

% create the exponent matrix that corresponds to the quadratic map
counter = 1;

for j = 1:N2
    Equad(:, counter:counter+N1-1) = Eext1 + Eext2(:, j) * ones(1, N1);
    counter = counter + N1;
end

% create the exponent matrix that corresponds to the cubic map
for j = 1:N3
    Ecub(:, (j - 1)*M+1:j*M) = Equad + Eext3(:, j) * ones(1, M);
end

% loop over all dimensions
for i = 1:length(ind)

    % loop over all quadratic matrices: \sum_k (pZ' T_k pZ) * pZ_k
    for k = 1:length(ind{i})

        % quadratic evaluation
        quadMat = Gext1' * T{i, ind{i}(k)} * Gext2;
        quadVec = reshape(quadMat, 1, []);


        % cubic generator matrix: loop over all generators
        for j = 1:N3
            Gcub(i, (j - 1)*M+1:j*M) = Gcub(i, (j - 1)*M+1:j*M) + ...
                quadVec * Gext3(ind{i}(k), j);
        end
    end
end

% add up all generators that belong to identical exponents
[ExpNew, Gnew] = removeRedundantExponents(Ecub, Gcub);

% mixed multiplication with the zonotope from the independent terms
zonoInd = cubMap(Z1, Zrem2, Z3, T, ind) + ...
    cubMap(Zrem1, Z2, Z3, T, ind) + ...
    cubMap(Zrem1, Zrem2, Z3, T, ind) + ...
    cubMap(Z1, Z2, Zrem3, T, ind) + ...
    cubMap(Z1, Zrem2, Zrem3, T, ind) + ...
    cubMap(Zrem1, Z2, Zrem3, T, ind) + ...
    cubMap(Zrem1, Zrem2, Zrem3, T, ind);

c = zonoInd.c;
GI = zonoInd.G;

% remove zero generators
GI(:, sum(abs(GI), 1) == 0) = [];

% construct the resulting polynomial zonotope
if all(ExpNew(:, 1) == 0)
    res = polyZonotope(c+Gnew(:, 1), Gnew(:, 2:end), GI, ...
        ExpNew(:, 2:end), id);
else
    res = polyZonotope(c, Gnew, GI, ExpNew, id);
end
end

% ------------------------------ END OF CODE ------------------------------
