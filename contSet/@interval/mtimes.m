function res = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for intervals
%
% Syntax:
%    res = factor1 * factor2
%    res = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - interval object, numeric
%    factor2 - interval object, numeric, contSet object
%
% Outputs:
%    res - interval
%
% Example:
%    factor1 = interval(-2,0);
%    factor2 = interval(1,2);
%    factor1 * factor2
%
% References:
%    [1] M. Althoff. "Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       19-June-2015
% Last update:   25-June-2015
%                18-November-2015
%                01-February-2016 (DG, Fixed a matrix case)
%                27-February-2016 (DG, New matrix case)
%                21-July-2016 (MA, case scalar and numeric added)
%                22-July-2016 (MA, case that factor1 is numeric has been added)
%                26-July-2016 (multiplication with zonotope added)
%                05-August-2016 (simplified some cases; matrix case corrected)
%                23-June-2022 (VG, added support for empty matrices/intervals)
%                04-April-2023 (TL, vectorized matrix case)
%                10-October-2023 (TL, fix for scalar 0*[inf,-inf])
%                04-October-2024 (TL, fix for matrix 0*[inf,-inf])
% Last revision: 04-October-2024 (MW, remove InferiorClasses)

% ------------------------------ BEGIN CODE -------------------------------

%%% interval x contSet

% in the cases below, factor1 must be an interval object (otherwise, we
% would not enter this function)
if isa(factor2,'polyZonotope')
    res = aux_mtimes_polyZonotope(factor1,factor2);
    return
elseif isa(factor2,'zonotope')
    res = aux_mtimes_zonotope(factor1,factor2);
    return
elseif isa(factor2,'conZonotope')
    res = aux_mtimes_conZonotope(factor1,factor2);
    return
elseif isa(factor2,'zonoBundle')
    res = aux_mtimes_zonoBundle(factor1,factor2);
    return
end

% other contSet cases not supported
if isa(factor2,'contSet') && ~isa(factor2,'interval')
    throw(CORAerror('CORA:noops',factor1,factor2));
end


%%% interval x interval
    
% scalar case
if isscalar(factor1) && isscalar(factor2)
    res = aux_mtimes_scalar(factor1,factor2);
    return
end

if isscalar(factor1)  % factor2 is a matrix
    res = aux_mtimes_scalar_matrix(factor1,factor2);
    return
end

if isscalar(factor2)  % factor1 is a matrix
    res = aux_mtimes_matrix_scalar(factor1,factor2);
    return
end

% matrix case
if ~issparse(factor1) && ~issparse(factor2)
    res = aux_mtimes_nonsparse(factor1,factor2);
    return
else
    res = aux_mtimes_sparse(factor1,factor2);
    return
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_mtimes_scalar(factor1,factor2)

% obtain possible values
if isnumeric(factor1)
    res = factor2;
    if factor1 == 0
        % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
        possibleValues = 0;
    else
        possibleValues = [factor1*factor2.inf, factor1*factor2.sup];
    end
elseif isnumeric(factor2)
    res = factor1;

    if factor2 == 0
        % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
        possibleValues = 0;
    else
        possibleValues = [factor1.inf*factor2, factor1.sup*factor2];
    end
else
    res = factor1;
    possibleValues = [factor1.inf*factor2.inf, factor1.inf*factor2.sup, ...
        factor1.sup*factor2.inf, factor1.sup*factor2.sup];
end

%infimum
res.inf = min(possibleValues);

%supremum
res.sup = max(possibleValues);

end

function res = aux_mtimes_scalar_matrix(factor1,factor2)

%obtain possible values
if isnumeric(factor1) 
    if factor1 < 0
        %infimum and supremum
        res = interval( ...
            factor1*factor2.sup, ...
            factor1*factor2.inf ...
        );
    elseif factor1 > 0
        %infimum and supremum
        res = interval( ...
            factor1*factor2.inf, ...
            factor1*factor2.sup ...
        );
    else % factor1 == 0
        % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
        res = interval(zeros(size(factor2.inf)));
    end
else
    res = factor1.*factor2;
end

end

function res = aux_mtimes_matrix_scalar(factor1,factor2)

%obtain possible values
if isnumeric(factor2)
    if factor2 < 0
        %infimum and supremum
        res = interval( ...
            factor2*factor1.sup, ...
            factor2*factor1.inf ...
        );
    elseif factor2 > 0
        %infimum and supremum
        res = interval( ...
            factor2*factor1.inf, ...
            factor2*factor1.sup ...
        );
    else % factor2 == 0
        % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
        res = interval(zeros(size(factor1.inf)));
    end
else
    res = factor1.*factor2;
end

end

function res = aux_mtimes_nonsparse(factor1,factor2)

% compute fast algorithm
% [m, k] * [k, n] = [m, n]
% -> [m, k, 1] .* [1, k, n] = [m, k, n]

if isnumeric(factor2) 
    [m, ~] = size(factor1);
    extSize = [1, size(factor2)];
    factor2 = reshape(factor2, extSize);
    res = factor1 .* factor2;

    % fix 0*Inf=NaN cases (see scalar case above)
    res.inf(isnan(res.inf)) = 0;
    res.sup(isnan(res.sup)) = 0;

    % [m,k,n] -> [m,n]
    res.inf = sum(res.inf, 2);
    res.sup = sum(res.sup, 2);
    res = reshape(res, m, []);
    
else
    [m, ~] = size(factor1);
    extSize = [1, size(factor2)];
    factor2.inf = reshape(factor2.inf, extSize);
    factor2.sup = reshape(factor2.sup, extSize);
    res = factor1 .* factor2;

    % fix 0*Inf=NaN cases (see scalar case above)
    res.inf(isnan(res.inf)) = 0;
    res.sup(isnan(res.sup)) = 0;

    % [m,k,n] -> [m,n]
    res.inf = sum(res.inf, 2);
    res.sup = sum(res.sup, 2);
    res = reshape(res, m, []);

end

end

function res = aux_mtimes_sparse(factor1,factor2)

% sparse only supports 2d arrays, keeping old algorithm..

% convert both to interval
factor1 = interval(factor1);
factor2 = interval(factor2);

I1 = factor1.inf;
S1 = factor1.sup;

[m, ~] = size(I1);
[~, n] = size(factor2.inf);

% preallocate output bounds [m,n]
resInf = zeros(m,n);
resSup = zeros(m,n);

A = interval(zeros(m,1));
for i = 1:m
    % get i-th row [1, k], transpose
    A.inf = I1(i, :)';
    A.sup = S1(i, :)';

    % [k, 1] .* [k, n] = [k, n]
    B = A .* factor2;
    resInf(i, :) = sum(B.inf, 1);
    resSup(i, :) = sum(B.sup, 1);
end
res = interval(resInf, resSup);

end

function Z = aux_mtimes_zonotope(I,Z)
% see Theorem 3.3 in [1]
% note: factor1 must be an interval object

% get center and radius of interval matrix
T = center(I);
S = rad(I);

Zabssum = sum(abs([Z.c,Z.G]),2);

% compute new zonotope
Z.c = T*Z.c;
Z.G = [T*Z.G,diag(S*Zabssum)];

end

function pZ = aux_mtimes_polyZonotope(I,pZ)
% note: factor1 must be an interval object

% get center and radius of interval
m = center(I);
r = rad(I);

% calculate interval over-approximation
I = interval(pZ);
s = abs(center(I)) + rad(I);

% compute new polyZonotope
pZ.c = m*pZ.c;
if ~isempty(pZ.G)
    pZ.G = m*pZ.G;
end
if ~isempty(pZ.GI)
    pZ.GI = [m*pZ.GI, diag(r*s)];
else
    pZ.GI = diag(r*s);
end
% pZ.id stays the same

end

function cZ = aux_mtimes_conZonotope(I,cZ)
    
% center and radius of interval matrix
m = center(I);
r = rad(I);

% absolute value of zonotope center and generators
Zabssum = sum(abs([cZ.c,cZ.G]),2);

% construct resulting conZonotope object
cZ.c = m*cZ.c;
cZ.G = [m*cZ.G, diag(r*Zabssum)];
cZ.A = [cZ.A, zeros(size(cZ.A,1),size(Zabssum,1))];
% cZ.b stays the same

end

function zB = aux_mtimes_zonoBundle(I,zB)

% loop over all individual zonotopes
for i=1:zB.parallelSets
    zB.Z{i} = aux_mtimes_zonotope(I,zB.Z{i});
end

end

% ------------------------------ END OF CODE ------------------------------
