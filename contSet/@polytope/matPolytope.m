function matP = matPolytope(P)
% matPolytope - converts the given polytope into a matPolytope object
%
% Syntax:
%    matP = matPolytope(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    matP - matPolytope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       21-June-2010
% Last update:   02-May-2024 (TL, moved out of constructor)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get vertices from polytope class
V=vertices(input);

% rewrite vertices as matrices (in reverse to pre-allocate entire matrix)
for i=length(V(:,1):-1:1)
    matrixVertex(:,:,i) = vec2mat(V(i,:));
end

matP = matPolytope(P);

end

% ------------------------------ END OF CODE ------------------------------
