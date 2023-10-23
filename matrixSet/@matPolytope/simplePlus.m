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

% Authors:       Matthias Althoff
% Written:       06-July-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize potential vertices
Vpot=[];
%Calculate posiible new vertices by adding all combinations
for j=1:summand1.verts
    for i=1:summand2.verts
        Vpot{end+1}=summand1.vertex{j}+summand2.vertex{i};
    end
end

%instantiate matrix polytope
matP=matPolytope(Vpot);

% ------------------------------ END OF CODE ------------------------------
