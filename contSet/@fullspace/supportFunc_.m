function [val,x] = supportFunc_(fs,dir,type,varargin)
% supportFunc_ - calculates the upper or lower bound of a full-dimensional
%    space along a certain direction
%    case R^0: 'upper' -> +Inf, 'lower' -> -Inf
%
% Syntax:
%    val = supportFunc_(fs,dir,type)
%
% Inputs:
%    fs - fullspace object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the full-dimensional space in the specified direction
%    x - support vector
%
% Example: 
%    fs = fullspace(2);
%    dir = [1;1];
%    [val,x] = supportFunc(fs,dir);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (rename supportFunc_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if fs.dimension == 0 && nargout == 2
    throw(CORAerror('CORA:notSupported',...
        'Intersection check of R^0 not supported'));
end

% set is always unbounded
if strcmp(type,'upper')
    val = Inf;
    if nargout == 2
        x = Inf(fs.dimension,1) .* sign(dir);
    end

elseif strcmp(type,'lower')
    val = -Inf;
    if nargout == 2
        x = -Inf(fs.dimension,1) .* sign(dir);
    end

elseif strcmp(type,'range')
    val = interval(-Inf,Inf);
    if nargout == 2
        x = [-Inf(fs.dimension,1),Inf(fs.dimension,1)] .* sign(dir);
    end

end

% ------------------------------ END OF CODE ------------------------------
