function [val,x] = supportFunc(hyp,dir,varargin)
% supportFunc - computes support function of constrained hyperplane in a 
%    given direction
%
% Syntax:  
%    [val,x] = supportFunc(hyp,dir)
%    [val,x] = supportFunc(hyp,dir,type)
%
% Inputs:
%    hyp - conHyperplane object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of H in the specified direction
%    x   - point for which holds: dir'*x=val
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    supportFunc(hyp,[-1;1]/sqrt(2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      22-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% convert to polytope and call support function for polytopes
if nargout <= 1
    val = supportFunc(mptPolytope(hyp),dir,varargin{:});
elseif nargout == 2
    [val,x] = supportFunc(mptPolytope(hyp),dir,varargin{:});
end

%------------- END OF CODE --------------