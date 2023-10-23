function [val,x] = supportFunc_(O,dir,type,varargin)
% supportFunc_ - calculates the upper or lower bound of an empty set along
%    a certain direction
%
% Syntax:
%    val = supportFunc_(O,dir,type)
%
% Inputs:
%    O - emptySet object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the full-dimensional space in the specified direction
%    x - support vector
%
% Example: 
%    O = emptySet(2);
%    dir = [1;1];
%    [val,x] = supportFunc(O,dir);
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

% bounds are fixed by theory
if strcmp(type,'upper')
    val = -Inf;
    x = double.empty(O.dimension,0);

elseif strcmp(type,'lower')
    val = Inf;
    x = double.empty(O.dimension,0);

elseif strcmp(type,'range')
    val = interval(-Inf,Inf);
    % actually, there would have to be two empty n-dimensional vectors next
    % to each other, but this is not supported...
    x = double.empty(O.dimension,0);
end

% ------------------------------ END OF CODE ------------------------------
