function P = polytope(hyp)
% polytope - Converts a constrained hyperplane to a polytope
%
% Syntax:
%    P = polytope(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    P - polytope object
%
% Example:
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    P = polytope(hyp);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/polytope

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% conversion
if isempty(hyp.C)
    A = [hyp.a;-hyp.a];
    b = [hyp.b;-hyp.b];       
else
    A = [hyp.a;-hyp.a;hyp.C];
    b = [hyp.b;-hyp.b;hyp.d];
end

% instantiate polytope
P = polytope(A,b);

% ------------------------------ END OF CODE ------------------------------
