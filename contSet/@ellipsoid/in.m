function res = in(obj1,obj2)
% in - checks whether obj2 is contained in obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%
% Inputs:
%    obj1 - ellipsoid object 
%    obj2 - contSet object 
%
% Outputs:
%    res - result of containment check
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2));
%    zono = zonotope([0 1 0;0 1 1]);
%
%    in(E1,E2)
%    in(E1,zono)
%
%    figure
%    hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%
%    figure
%    hold on
%    plot(E1,[1,2],'b');
%    plot(zono,[1,2],'r');
%
% References:
%            [1] Yildirim, E.A., 2006. On the minimum volume covering 
%                ellipsoid of ellipsoids. SIAM Journal on Optimization, 
%                17(3), pp.621-641.     
%            [2] SDPT3: url: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann, Niklas Kochdumper
% Written:      15-October-2019 
% Last update:  21-November-2019 (NK, extend to other sets)
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;

    if isnumeric(obj2)
        
        for i = 1:size(obj2,2)
            res = containsPoint(obj1,obj2(:,i)); 
            if res ~= 1
                return;
            end
        end
        
    elseif isa(obj2,'ellipsoid')
        
        res = ellipsoidInEllipsoid(obj1,obj2);
        
    else
        
        % zonotope over-approximation
        if isa(obj2,'capsule') || isa(obj2,'taylm') || isa(obj2,'polyZonotope')
           obj2 = zonotope(obj2); 
        end
        
        % check if all vertices of the set are contained
        res = in(obj1,vertices(obj2));
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = containsPoint(E,p)
% check if a point is inside an ellipsoid

    p_ = p-E.q;
    y = linsolve(E.Q,p_);
    res = p_'*y <= 1;

end

function res = ellipsoidInEllipsoid(E1,E2)
% checks whether ellipsoid E2 is contained in ellipsoid E2

    %For more detail, see [1]
    t = sdpvar;
    condTOL = 1e6;
    if E2.dim ~=E1.dim
        error('dimensions do not match');
    elseif E2.dim<length(E2.Q) || E1.dim<length(E1.Q)
        error('E1 or E2 not full-dimensional');
    elseif cond(E2.Q)>condTOL || cond(E1.Q)>condTOL
        warning('One of the ellipsoids is very squished in at least one direction - computation might fail');
    end
    Q1 = E2.Q;
    Q2 = E1.Q;
    q1 = E2.q;
    q2 = E1.q;
    %follows from the fact that Q2-Q1>=0 (psd) together with Q1,Q2 invertible
    %implies inv(A)>=inv(B).
    if q1==q2
        res = min(eig(Q2-Q1))>=0;
        return;
    end
    if ~isYalmipInstalled()
        error('YALMIP must be on the MATLAB search path to use this function');
    end
    Q1inv = inv(Q1);
    Q2inv = inv(Q2);
    M1 = [Q1inv, -Q1inv*q1;
          -(Q1inv*q1)', q1'*Q1inv*q1-1];
    M2 = [Q2inv, -Q2inv*q2;
          -(Q2inv*q2)', q2'*Q2inv*q2-1];

    options = sdpsettings('verbose',0);
    %when solving problems with higher-dimensional ellipsoids, sdpt3 [2] is
    %recommended to avoid numerical issues
    diagnostics = optimize(t*M1>=M2,[],options);
    %either feasible or not feasible
    if diagnostics.problem~=1 && diagnostics.problem~=0
        error('error with solving feasibility problem!');
    end
    if value(t)==0
        error('problem with solution of feasibility problem (t)');
    end
    res = ~diagnostics.problem;
end

%------------- END OF CODE --------------