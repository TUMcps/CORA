function val = supportFunc_(tay,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of a Taylor model along 
%    a certain direction
%
% Syntax:
%    val = supportFunc_(tay,dir)
%    val = supportFunc_(tay,dir,type)
%
% Inputs:
%    tay - taylm object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the Taylor model in the specified direction
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, conZonotope/supportFunc_

% Authors:       Niklas Kochdumper
% Written:       19-November-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------
    
% compute enclosing interval of Taylor model projeted onto the direction
I = interval(dir'*tay);

% upper or lower bound
if strcmp(type,'lower')
    val = infimum(I);
elseif strcmp(type,'range')
    val = supremum(I);
else
    % 'val' is already the desired result
end

% ------------------------------ END OF CODE ------------------------------
