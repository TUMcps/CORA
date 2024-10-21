function S_out = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or a 
%    matrix zonotope with a matrix zonotope
%
% Syntax:
%    matZ = factor1 * factor2
%    matZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numeric matrix or matZonotope object
%    factor2 - numeric matrix or matZonotope object
%
% Outputs:
%    matZ - matrix zonotope
%
% Example: 
%
% Reference:
%    [1] M. Althoff et al. "Modeling, Design, and Simulation of Systems 
%        with Uncertainties". 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       18-June-2010 
% Last update:   05-August-2010
%                25-April-2024 (TL, pagemtimes, much faster implementation)
% Last revision: 04-October-2024 (MW, remove InferiorClasses from contSet)

% ------------------------------ BEGIN CODE -------------------------------

% factor1 is a numeric matrix -> factor2 must be a matZonotope object
if isnumeric(factor1)
    % copy matZonotope object
    S_out=factor2;
    % evaluate linear map
    S_out.C = factor1*factor2.C;
    S_out.G = pagemtimes(factor1,factor2.G);
    return
end

% factor2 is a numeric matrix -> factor1 must be a matZonotope object
if isnumeric(factor2)
    % copy intervalMatrix object
    S_out = factor1;
    % evaluate linear map
    S_out.C = factor1.C*factor2;
    S_out.G = pagemtimes(factor1.G,factor2);
    return
end

if isa(factor2,'matZonotope')
    S_out = aux_mtimes_matZonotope(factor1,factor2);
    return
end


%%% below, multiplications of matZonotope * contSet

if isa(factor2,'zonotope')
    S_out = aux_mtimes_zonotope(factor1,factor2);
    return
end

if isa(factor2,'zonoBundle')
    S_out = aux_mtimes_zonoBundle(factor1,factor2);
    return
end

if isa(factor2,'polyZonotope')
    % not implemented yet...
    % see: Ertai Luo, Niklas Kochdumper, and Stanley Bak. "Reachability
    % Analysis for Linear Systems with Uncertain Parameters using Polynomial
    % Zonotopes," In Proceedings of the ACM International Conference on
    % Hybrid Systems: Computation and Control. Article 17, 1–12, 2023.
end

throw(CORAerror('CORA:noops',factor1,factor2));

end


% Auxiliary functions -----------------------------------------------------

function S_out = aux_mtimes_matZonotope(factor1,factor2)

% initialize matrix zonotope 1
matZ1=factor1;
% initialize matrix zonotope 2
matZ2=factor2;

% concat center with generators
Z1 = cat(3,matZ1.C,matZ1.G);
[n1,m1,h1] = size(Z1);
Z2 = cat(3,matZ2.C,matZ2.G);
[n2,m2,h2] = size(Z2);

% reshape Z2
Z2 = reshape(Z2,n2,m2,1,h2);

% multiply each matrix from either set
Z = pagemtimes(Z1,Z2);

% reshape back ---

% fix dimensions for scalar multiplication
if n1 == 1 && m1 == 1
    n1 = n2;
end
if n2 == 1 && m2 == 1
    m2 = m1;
end

%  reshape
Z = reshape(Z,n1,m2,h1*h2);

% construct final matrix zonotope
S_out = matZonotope();
S_out.C = Z(:,:,1);
S_out.G = Z(:,:,2:end);

end

function Z = aux_mtimes_zonotope(matZ,Z)
% see Sec. 4.4.1 in [1]

% extract center and generators
c = Z.c;
G = Z.G;

% get output dimension
if all(dim(matZ) == 1)
    % multiplication with a scalar
    n = length(c);
else
    % matrix multiplication
    n = matZ.dim(1);
end

% obtain first zonotope
Z.c = matZ.C * c;
Z.G = [
    matZ.C * G, ...
    reshape(pagemtimes(matZ.G,c),n,[]), ...
    reshape( ...
        pagemtimes(matZ.G,reshape(G,size(G,1),1,1,size(G,2))), ...
        n,[])
];

end

function zB = aux_mtimes_zonoBundle(matZ,zB)
% compute linear map for each zonotope

for i=1:zB.parallelSets
    zB.Z{i} = aux_mtimes_zonotope(matZ,zB.Z{i});
end

end

% ------------------------------ END OF CODE ------------------------------
 