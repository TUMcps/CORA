function S_out = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of matrix or an 
%    interval matrix with an interval matrix
%
% Syntax:
%    intMat = factor1 * factor2
%    intMat = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - intervalMatrix object, numerical matrix
%    factor2 - intervalMatrix object, numerical matrix, contSet object
%
% Outputs:
%    S_out - result of the linear map (intervalMatrix or contSet object)
%
% Example:
%
% References:
%    [1] M. Althoff. "Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/mtimes

% Authors:       Matthias Althoff
% Written:       18-June-2010 
% Last update:   05-August-2010
%                03-April-2022 (MW, remove setting check)
% Last revision: 04-October-2024 (MW, remove InferiorClasses from contSet)

% ------------------------------ BEGIN CODE -------------------------------

% factor1 is a numeric matrix -> factor2 must be an intervalMatrix object
if isnumeric(factor1)
    % copy intervalMatrix object
    S_out = factor2;
    % evaluate linear map
    S_out.int = factor1*factor2.int;
    return
end
    
% factor2 is a numeric matrix -> factor1 must be an intervalMatrix object
if isnumeric(factor2)
    % copy intervalMatrix object
    S_out = factor1;
    % evaluate linear map
    S_out.int = factor1.int*factor2;
    return
end

if isa(factor2,'intervalMatrix')
    % copy intervalMatrix object
    S_out = factor1;
    % evaluate linear map
    S_out.int = factor1.int*factor2.int;
    return
end


%%% below, multiplications of intervalMatrix * contSet

if isa(factor2,'zonotope')
    S_out = aux_mtimes_zonotope(factor1,factor2);
    return
end

if isa(factor2,'polyZonotope')
    S_out = aux_mtimes_polyZonotope(factor1,factor2);
    return
end

if isa(factor2,'zonoBundle')
    S_out = aux_mtimes_zonoBundle(factor1,factor2);
    return
end

if isa(factor2,'emptySet')
    S_out = emptySet(dim(factor1,1));
    return
end

throw(CORAerror('CORA:noops',factor1,factor2));

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_mtimes_zonotope(intMat,Z)
% see Theorem 3.3 in [1]

% get minimum and maximum
M_min = infimum(intMat.int);
M_max = supremum(intMat.int); 

% get center of interval matrix
T = 0.5*(M_max+M_min);
% get symmetric interval matrix
S = 0.5*(M_max-M_min);
Zabssum = sum(abs([Z.c,Z.G]),2);

Z.c = T*Z.c;
Z.G = [T*Z.G,diag(S*Zabssum)];

end

function pZ = aux_mtimes_polyZonotope(intMat,pZ)

% get minimum and maximum
M_min = infimum(intMat.int);
M_max = supremum(intMat.int); 

% get center and radius of interval matrix
M = 0.5*(M_max+M_min);
R = 0.5*(M_max-M_min);

% calculate interval over-approximation
I = interval(pZ);
S = abs(center(I)) + rad(I);

% compute new polyZonotope
pZ.c = M*pZ.c;
if ~isempty(pZ.G)
    pZ.G = M*pZ.G;
end
if ~isempty(pZ.GI)
    pZ.GI = [M*pZ.GI, diag(R*S)];
else
    pZ.GI = diag(R*S);
end
% pZ.id stays the same

end

function zB = aux_mtimes_zonoBundle(intMat,zB)
% compute linear map for each zonotope

for i=1:zB.parallelSets
    zB.Z{i} = aux_mtimes_zonotope(intMat,zB.Z{i});
end

end

% ------------------------------ END OF CODE ------------------------------
