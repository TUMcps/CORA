function [val,x] = supportFunc_(I,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of an interval along a
%    certain direction
%
% Syntax:
%    val = supportFunc_(I,dir,type)
%    [val,x] = supportFunc_(I,dir,type)
%
% Inputs:
%    I - interval object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the interval in the specified direction
%    x - support vector
%
% Example:
%    I = interval([-2;1],[3;2]);
%    dir = [1;1]/sqrt(2);
%    supportFunc(I,dir)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, zonotope/supportFunc_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       19-November-2019
% Last update:   27-March-2023 (MW, rename supportFunc_)
% Last revision: 06-April-2023 (MW, rewrite function)

% ------------------------------ BEGIN CODE -------------------------------

% special handling for empty set
if representsa_(I,'emptySet',0)
    x = [];
    if strcmp(type,'upper')
        val = -Inf;
    elseif strcmp(type,'lower')
        val = Inf;
    elseif strcmp(type,'range')
        val = interval(-Inf,Inf);
    end
    return
end

% take infimum/supremum depending on sign of direction; for entries with 0,
% it does not matter
idx = sign(dir) == -1;
if strcmp(type,'upper')
    x = I.sup;
    x(idx) = I.inf(idx);
    val = dir'*x;
elseif strcmp(type,'lower')
    x = I.inf;
    x(idx) = I.sup(idx);
    val = dir'*x;
elseif strcmp(type,'range')
    x = [I.inf I.sup];
    x(idx,1) = I.sup(idx);
    x(idx,2) = I.inf(idx);
    val = interval(dir'*x(:,1),dir'*x(:,2));
end

% % old version:
% if nargout == 1
%     val = supportFunc_(zonotope(I),dir,varargin{:});
% else
%     [val,x] = supportFunc_(zonotope(I),dir,varargin{:});
% end

% ------------------------------ END OF CODE ------------------------------
