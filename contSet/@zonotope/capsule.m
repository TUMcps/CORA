function C = capsule(Z)
% capsule - encloses a zonotope with a capsule
%
% Syntax:
%    C = capsule(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    C - capsule object
%
% Example: 
%    Z = zonotope([1;-1],[2 -3 1; 0.5 1 -2]);
%    C = capsule(Z);
%
%    figure; hold on;
%    plot(Z,[1,2],'r');
%    plot(C,[1,2],'b');
%
% Other m-files required: interval (constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       23-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute orthogonal basis using PCA
G = generators(Z);
[B,~,~] = svd([-G,G]);

% compute enclosing interval in the transformed space
int = interval(B'*Z);

% enclose interval with a capsule
C = capsule(int);

% back-transformation to original space
C = B*C;

% ------------------------------ END OF CODE ------------------------------
