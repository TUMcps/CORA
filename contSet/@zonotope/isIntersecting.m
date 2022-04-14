function res = isIntersecting(Z,obj,varargin)
% isIntersecting - determines if zonotope Z intersects obj
%
% Syntax:  
%    res = isIntersecting(Z,obj)
%    res = isIntersecting(Z,obj,type)
%
% Inputs:
%    Z - zonotope object
%    obj - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    zono1 = zonotope([0 1 1 0;0 1 0 1]);
%    zono2 = zonotope([2 -1 1 0;2 1 0 1]);
%    zono3 = zonotope([3.5 -1 1 0;3 1 0 1]);
%
%    isIntersecting(zono1,zono2)
%    isIntersecting(zono1,zono3)
%
%    figure; hold on;
%    plot(zono1,[1,2],'b');
%    plot(zono2,[1,2],'g');
%
%    figure; hold on;
%    plot(zono1,[1,2],'b');
%    plot(zono3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isIntersecting

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % get zonotope object
    if ~isa(Z,'zonotope')
        % switch variables
        temp = Z;
        Z = obj;
        obj = temp;
    end
    
    % call function for other set representations
    if isa(obj,'halfspace') || isa(obj,'conHyperplane') || ...
       isa(obj,'mptPolytope') || isa(obj,'ellipsoid')
   
        res = isIntersecting(obj,Z,varargin{:});
   
    else
        
        res = isIntersecting(conZonotope(Z),obj,varargin{:});
        
    end     
end

%------------- END OF CODE --------------