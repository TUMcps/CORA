function [val,x,fac] = supportFunc_(pgon,dir,type,varargin)
% supportFunc_ - calculates the upper or lower bound of a polygon along a
%    certain direction
%
% Syntax:
%    [val,x,fac] = supportFunc_(pgon,dir,type)
%
% Inputs:
%    pgon - polygon object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the polygon in the specified direction
%    x - support vector
%    fac - factor values that correspond to the upper bound
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init values
x = []; fac = [];

% return correct values for empty polygon
if representsa_(pgon,"emptySet",eps)
    if strcmp(type,'upper')
        val = -Inf;
    elseif strcmp(type,'lower')
        val = Inf;
    elseif strcmp(type,'range')
        val = interval(-Inf,Inf);
    end
    return
end

% get vertices
V = vertices_(pgon);

% compute support func
s = dir' * V;

% upper or lower bound
if strcmp(type,'lower')
    [val,idx] = min(s);
    x = V(:,idx);
    
elseif strcmp(type,'upper')
    [val,idx] = max(s);
    x = V(:,idx);
    
elseif strcmp(type,'range')
    % compute lower/upper
    [val_lower,x_lower] = supportFunc_(pgon,dir,'lower');
    [val_upper,x_upper] = supportFunc_(pgon,dir,'upper');

    % obtain result
    val = interval(val_lower,val_upper);
    x = [x_lower x_upper];
end


% ------------------------------ END OF CODE ------------------------------
