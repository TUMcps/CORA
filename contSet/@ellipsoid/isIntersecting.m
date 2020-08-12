function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if ellipsoid obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - ellipsoid object
%    obj2 - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2),[1;0]);
%    E3 = ellipsoid(0.3*eye(2),[2;0]);
%
%    isIntersecting(E1,E2)
%    isIntersecting(E1,E3)
%
%    figure
%    hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%
%    figure
%    hold on
%    plot(E1,[1,2],'b');
%    plot(E3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isIntersecting

% Author:       Victor Gassmann, Niklas Kochdumper
% Written:      13-March-2019 
% Last update:  21-November-2019 (NK: extendet to other sets)
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end

    % get ellipsoid object
    if ~isa(obj1,'ellipsoid')
        temp = obj1;
        obj1 = obj2;
        obj2 = temp;
    end
    
    % check for intersection
    if isa(obj2,'ellipsoid')
       res = intersectionEllipsoid(obj1,obj2);
    elseif isa(obj2,'hyperplane')
       res = isIntersecting(obj2,obj1,type); 
    else
       
       % exact or over-approximative algorithm
       if strcmp(type,'exact')
           error('No exact algorithm implemented for this set representation!');
       else

           % zonotope over-approximation
           m = 5*obj1.dim-1;
           zono = zonotope(obj1,m,'o:norm');
           
           % check for intersection
           res = isIntersecting(zono,obj2); 
       end     
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = intersectionEllipsoid(E1,E2)

    %We assume that ellipsoids are full-dimensional
    if E1.isdegenerate || E2.isdegenerate
        error('E1/E2 degenerate, not supported');
    elseif length(E1.Q)~=length(E2.Q)
        error('E1/E2 have to have same dimensions');
    end
    %check if E1==E2
    if E1==E2
        res = true;
        return;
    end
    if ~isYalmipInstalled()
        error('YALMIP must be on the MATLAB search path to use this function');
    end
    Q1inv = inv(E1.Q);
    Q2inv = inv(E2.Q);
    q1 = E1.q;
    q2 = E2.q;
    x = sdpvar(length(E1.Q),1);
    %check if below constraint set is feasible => means E1,E2 intersecting
    F = [(x-q1)'*Q1inv*(x-q1)<=1, (x-q2)'*Q2inv*(x-q2)<=1];
    opts = sdpsettings('verbose',0);
    ds = optimize(F,[],opts);
    if ds.problem==0
        res = true;
    elseif any(ds.problem==[1,12,15])
        res = false;
    else
        error('error with solving feasibility problem!');
    end
end

%------------- END OF CODE --------------