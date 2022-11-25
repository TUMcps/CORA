function val = supportFunc(obj,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a taylm object
%               along a certain direction
%
% Syntax:  
%    val = supportFunc(obj,dir)
%    val = supportFunc(obj,dir,type)
%
% Inputs:
%    obj - taylm object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % parse input arguments
    type = 'upper';
    
    if nargin >= 3 && ~isempty(varargin{1})
        type = varargin{1};
    end

    % project Taylor model onto the direction
    obj_ = dir'*obj;
    
    % compute enclosing interval
    int = interval(obj_);
    
    % upper or lower bound
    if strcmp(type,'lower')
       val = infimum(int);
    else        
       val = supremum(int);
    end
end

%------------- END OF CODE --------------