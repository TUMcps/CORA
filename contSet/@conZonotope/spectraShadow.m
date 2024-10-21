function SpS = spectraShadow(cZ)
% spectraShadow - Converts a constrained zonotope to a spectrahedral shadow
%
% Syntax:
%    SpS = spectraShadow(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    SpS - spectraShadow object
%
% Example:
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%    SpS = spectraShadow(cZ);
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

G = cZ.generators;

if isempty(G)
    SpS = spectraShadow(1,cZ.c,zeros([size(cZ.c,1) 0]));
    return
end

m = size(G,2);

% We begin by setting up the block of constraints necessary to represent
% the hypercube [-1,1]^m
A0 = speye(2*m);
Ai = cell([1 m]);
for i = 1:m
    temp = zeros([2*m 1]);
    temp(i) = 1;
    temp(m+i) = -1;
    Ai{i} = spdiags(temp,0,2*m,2*m);
end

% Append constraints of the inequality constraints
if isempty(cZ.A)
    % If there are none, we are done      
else
    for j=1:size(cZ.A,1)
        A0 = blkdiag(A0, sparse([cZ.b 0; 0 -cZ.b]));
        for i=1:m
            Ai{i} = blkdiag(Ai{i}, sparse([-cZ.A(j,i) 0; 0 cZ.A(j,i)]));
        end
    end
end

% Concatenate everything
A = [A0 cat(2,Ai{:})];


% instantiate spectraShadow
SpS = spectraShadow(A, cZ.c, G);

% Additional properties
SpS.bounded.val = true;
SpS.emptySet.val = representsa_(cZ,'emptySet',1e-10);
SpS.fullDim.val = isFullDim(cZ);
SpS.center.val = center(cZ);

% ------------------------------ END OF CODE ------------------------------
