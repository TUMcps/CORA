function [val,x] = supportFunc_(hyp,dir,varargin)
% supportFunc_ - computes support function of constrained hyperplane in a 
%    given direction
%
% Syntax:
%    [val,x] = supportFunc_(hyp,dir)
%    [val,x] = supportFunc_(hyp,dir,type)
%
% Inputs:
%    hyp - conHyperplane object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the hyperplane in the specified direction or interval
%          containing both bounds
%    x - point(s) for which holds: dir'*x=val
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    supportFunc(hyp,[-1;1]/sqrt(2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc

% Authors:       Victor Gassmann
% Written:       22-March-2021
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------

% convert to polytope and call support function for polytopes
if nargout <= 1
    val = supportFunc_(polytope(hyp),dir,varargin{:});
elseif nargout == 2
    [val,x] = supportFunc_(polytope(hyp),dir,varargin{:});
end

% ------------------------------ END OF CODE ------------------------------
