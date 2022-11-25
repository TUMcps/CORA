function [val,x] = supportFunc(obj,d,type)
% supportFunc - computes support function of constrained hyperplane in
% direction d
%
% Syntax:  
%    [val,x] = supportFunc(obj,d)
%    [val,x] = supportFunc(obj,d,type)
%
% Inputs:
%    obj - conHyperplane object
%    d   - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of H in the specified direction
%    x   - point for which holds: dir'*x=val
%
% Example: 
%    ---
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
if ~exist('type','var')
    type = 'upper';
end

if ~strcmp(type,'upper') && ~strcmp(type,'lower')
    error('Third argument can either be "upper" or "lower"!');
end

if ~isa(d,'double') || ~any(size(d)==1) || length(d)~=dim(obj)
    error('Second argument has to be vector of appropriate dimension!');
end

obj = mptPolytope(obj);
if nargout<=1
    val = supportFunc(obj,d,type);
elseif nargout==2
    [val,x] = supportFunc(obj,d,type);
end
%------------- END OF CODE --------------