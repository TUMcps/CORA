function matP = simplePlus(summand1,summand2)
% simplePlus - computes the Minkowski addition of two matrix polytopes 
%    without reducing the vertices by a convex hull computation
%
% Syntax:
%    matP = simplePlus(summand1,summand2)
%
% Inputs:
%    summand1 - matPolytope object
%    summand2 - matPolytope object
%
% Outputs:
%    matP - matrix polytope after Minkowski addition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       06-July-2010 
% Last update:   02-May-2024 (TL, new structure of V)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

  % initialize matrix vertices
matV1 = summand1.V;
[n1,m1,h1] = size(Z1);
matV2 = summand2.V;
[n2,m2,h2] = size(matV2);

% reshape V2
matV2 = reshape(matV2,n2,m2,1,h2);

% obtain all potential vertices by adding all combinations
matV = matV1 + matV2;

% rewrite result as a matrix polytope
matV = reshape(matV,max(n1,n2),max(m1,m2),h1*h2);

% init matPolyzonotope
matP=matPolytope(matV);

end

% ------------------------------ END OF CODE ------------------------------
