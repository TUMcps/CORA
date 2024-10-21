function SpS = spectraShadow(P)
% spectraShadow - Converts a polytope to a spectrahedral shadow
%
% Syntax:
%    SpS = spectraShadow(P)
%
% Inputs:
%    P - P object
%
% Outputs:
%    SpS - spectraShadow object
%
% Example:
%    A = [1 2; -1 1; -1 -3; 2 -1];
%    b = ones(4,1);
%    P = polytope(A,b);
%    SpS = spectraShadow(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Deduce dimension of ambient space
n = dim(P);

% Define the inequality constraints Cx <= d that define P; this includes 2
% inequality constraints per equality constraint (also, exceptionally we
% use the notation (C,d) instead of (A,b) in order to avoid conflicting
% with the notation for spectrahedra :)

% compute constraints if not available
constraints(P);

% read constraints
C = [P.A;P.Ae;-P.Ae];
d = [P.b;P.be;-P.be];

% Initialize constraints of the equality constraint
A0 = [];
Ai = cell([1 n]);

% Append constraints of the inequality constraints
if isempty(C)
    % If there are none, we need to manually define a fullspace
    % spectrahedral shadow
    SpS = spectraShadow.Inf(n);
else
    for j=1:size(C,1)
        A0 = blkdiag(A0, [d(j)]);
        for i=1:n
            Ai{i} = blkdiag(Ai{i}, [-C(j,i)]);
        end
    end
    
    % Concatenate everything
    A = [A0 cat(2,Ai{:})];

    % instantiate spectraShadow
    SpS = spectraShadow(A);
    
    % Additional properties
    SpS.bounded.val = isBounded(P);
    SpS.emptySet.val = representsa_(P,'emptySet',1e-10);
    SpS.fullDim.val = isFullDim(P);
    SpS.center.val = center(P);
end

% ------------------------------ END OF CODE ------------------------------
