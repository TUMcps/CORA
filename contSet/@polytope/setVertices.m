function setVertices(P,V)
% setVertices - adds vertices to a polytope; caution: the correctness of
%    the vertices w.r.t the constraints has to be ensured by the user
%
% Syntax:
%    setVertices(P,V)
%
% Inputs:
%    P - polytope object
%    V - vertices
%
% Outputs:
%    P - polytope object with V in P.V
%
% Example: 
%    P = polytope([1 0; -1 1; -1 -1],[1;1;1]);
%    setVertices(P,[-1 0; 1 -2; 1 2]');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       13-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% insert vertices into the representation
P.V_.val = V;
P.isVRep.val = true;

% ------------------------------ END OF CODE ------------------------------
