function res = eq(hyp1,hyp2,varargin)
% eq - overloaded '==' operator for exact comparison of two hyperplanes
%
% Syntax:
%    res = hyp1 == hyp2
%    res = eq(hyp1,hyp2)
%    res = eq(hyp1,hyp2,tol)
%
% Inputs:
%    hyp1 - hyperplane object
%    hyp2 - hyperplane object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp1 = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    hyp2 = conHyperplane(halfspace([1;-1],0),[1 0;-1 0],[2;2]);
%    hyp1 == hyp2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/isequal

% Authors:       Mingrui Wang
% Written:       21-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(hyp1,hyp2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
