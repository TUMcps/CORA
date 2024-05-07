function matP = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of two matrix
%    polytopes or a matrix polytope with a matrix
%
% Syntax:
%    matP = plus(summand1,summand2)
%
% Inputs:
%    summand1 - matPolytope object or numerical matrix
%    summand2 - matPolytope object or numerical matrix
%
% Outputs:
%    matP - matrix polytope after Minkowski addition
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       21-June-2010 
% Last update:   02-May-2024 (TL, new structure of V)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%f ind a matrix polytope object
[matP,summand] = findClassArg(summand1,summand2,'matPolytope');

% is summand a matrix polytope?
if isa(summand,'matPolytope')

    % initialize matrix vertices
    matV1 = matP.V;
    [n1,m1,h1] = size(matV1);
    matV2 = summand.V;
    [n2,m2,h2] = size(matV2);

    % reshape V2
    matV2 = reshape(matV2,n2,m2,1,h2);

    % obtain all potential vertices by adding all combinations
    matVpot = matV1 + matV2;

    % reshape to vector for convex hull
    vecVpot = reshape(matVpot,max(n1,n2)*max(m1,m2),h1*h2);
   
    try
        % compute convex hull
        opt{1}='QJ';
        K=convhulln(vecVpot, opt);

        % only choose points on boundary
        vecV = vecVpot(K,:);
    catch
        CORAwarning('CORA:matrixSets','Convex hull could not be computed. Continuing with all vertices.')
        vecV = vecVpot;
    end

    % rewrite result as a matrix polytope
    matV = reshape(vecV,max(n1,n2),max(m1,m2),h1*h2);

    % init matPolyzonotope
    matP=matPolytope(matV);
    
% is summand a matrix?
elseif isnumeric(summand)
    % calculate minkowski sum
    matP.V = matP.V + summand;
end

% ------------------------------ END OF CODE ------------------------------
