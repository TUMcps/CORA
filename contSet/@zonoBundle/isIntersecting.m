function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if zonotope bundle obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - zonoBundle object
%    obj2 - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    int1 = interval([2;2],[4;4]);
%    int2 = interval([3.5;3],[5;5]);
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    isIntersecting(zB,int1)
%    isIntersecting(zB,int2)
%
%    figure; hold on
%    plot(zB,[1,2],'b');
%    plot(int1,[1,2],'g');
%    
%    figure; hold on
%    plot(zB,[1,2],'b');
%    plot(int2,[1,2],'r');
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

    % get zonoBundle object
    if ~isa(obj1,'zonoBundle')
        temp = obj1;
        obj1 = obj2;
        obj2 = temp;
    end
    
    % call function for other set representations
    if isa(obj2,'halfspace') || isa(obj2,'conHyperplane') || ...
       isa(obj2,'mptPolytope') || isa(obj2,'ellipsoid')
   
        res = isIntersecting(obj2,obj1,varargin{:});
   
    else
        
        res = isIntersecting(conZonotope(obj1),obj2,varargin{:});
        
    end     
end

%------------- END OF CODE --------------