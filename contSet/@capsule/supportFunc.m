function [val,x] = supportFunc(obj,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a capsule object
%            along a certain direction
%
% Syntax:  
%    [val,x] = supportFunc(obj,dir)
%    [val,x] = supportFunc(obj,dir,type)
%
% Inputs:
%    obj - capsule object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%    x - support vector
%
% Example: 
%    C = capsule([1; 1], [0.5; -1], 0.5);
%
%    val = supportFunc(C,[1;1]);
%   
%    figure
%    hold on
%    plot(C,[1,2],'r');
%    plot(conHyperplane(halfspace([1;1],val),[],[]),[1,2],'g');
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
    c = obj.c;
    g = obj.g;
    r = obj.r;
    
    % consider case where capsule is just a ball
    if isempty(g)
       g = zeros(size(c)); 
    end

    % compute upper or lower bound
    if strcmp(type,'upper')
       val = dir'*c + abs(dir'*g) + r*norm(dir);  
       x = c + g*sign(dir'*g) + r*dir;
    else
       val = dir'*c - abs(dir'*g) - r*norm(dir);  
       x = obj.c - g*sign(dir'*g) - r*dir;
    end
end

%------------- END OF CODE --------------