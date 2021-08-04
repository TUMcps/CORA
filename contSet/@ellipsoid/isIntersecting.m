function res = isIntersecting(obj,S,varargin)
% isIntersecting - determines if ellipsoid obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj,S)
%    res = isIntersecting(obj,S,type)
%
% Inputs:
%    obj - ellipsoid object
%    S - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2),[1;0]);
%    E3 = ellipsoid(0.3*eye(2),[2;0]);
%    Z = zonotope([3;0],0.5*eye(2));
%
%    isIntersecting(E1,E2)
%    isIntersecting(E1,E3)
%    isIntersecting(E1,Z,'approx')
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(E3,[1,2],'r');
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(Z,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isIntersecting

% Author:       Victor Gassmann, Niklas Kochdumper
% Written:      13-March-2019 
% Last update:  21-November-2019 (NK: extended to other sets)
%               10-March-2021 (refactored, simplified)
% Last revision:---

%------------- BEGIN CODE --------------
    if isempty(S)
        res = true;
        return;
    elseif isempty(obj)
        res = false;
        return;
    end
    
    % check if obj.Q is all-zero
    if rank(obj)==0
        res = in(S,obj.q);
        return;
    end
    
    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    if isa(S,'double')
        res = in(obj,S);
        return;
    end
    
    try
        % use distance
        res = distance(obj,S) <= obj.TOL;
    catch
        % distance not implemented for S
        % use fallback
        % exact or over-approximative algorithm
       if strcmp(type,'exact')
           error('No exact algorithm implemented for this set representation!');
       end
       % it is not certain that S implements quadMap
       if ismethod(S,'quadMap') && ismethod(obj,'interval')
           res = isIntersectingMixed(obj,S);
       else
           error('Not implemented for this set representation!');
       end
    end
end