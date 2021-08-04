function [val,x,fac] = supportFunc(obj,dir,varargin)
% supportFunc - calculates the upper or lower bound of a zonotope object
%               along a certain direction
%
% Syntax:  
%    val = supportFunc(obj,dir)
%    val = supportFunc(obj,dir,type)
%
% Inputs:
%    obj - zonotope object
%    dir - direction for which the bounds are calculated (vector)
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%    x - support vector
%    fac - factor values that correspond to the upper bound
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
    
    if size(dir,1) == 1 && size(dir,2) >= 1
        % transpose dir
        dir = dir';
    end
    
    if nargin >= 3 && ~isempty(varargin{1})
        type = varargin{1};
    end

    % project zonotope onto the direction
    obj_ = dir'*obj;

    % get object properties
    c = center(obj_);
    G = generators(obj_);
    
    % upper or lower bound
    if strcmp(type,'lower')
        
       val = c - sum(abs(G));
       fac = -sign(G)';
        
    else
        
       val = c + sum(abs(G));
       fac = sign(G)';
       
    end
    
    % compute support vector
    x = center(obj) + generators(obj)*fac;
end

%------------- END OF CODE --------------