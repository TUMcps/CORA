function [val,x] = supportFunc(obj,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a mptPolytope object
%            along a certain direction
%
% Syntax:  
%    val = supportFunc(obj,dir)
%    val = supportFunc(obj,dir,type)
%
% Inputs:
%    obj - mptPolytope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%    x - support vector
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

    % get object properties
    A = obj.P.A;
    b = obj.P.b;
    
    % linear program options
    options = optimoptions('linprog','display','off');
    
    % upper or lower bound
    if strcmp(type,'lower')
        
       % solve linear program
       [x,val] = linprog(dir',A,b,[],[],[],[],options);
        
    else
        
       % solve linear program
       [x,val] = linprog(-dir',A,b,[],[],[],[],options);
       val = -val;
       
    end
end

%------------- END OF CODE --------------