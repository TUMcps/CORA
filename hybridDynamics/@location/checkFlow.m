function [res,R] = checkFlow(obj,guard,R,options)
% checkFlow - remove all intersections for which the flow of the system
%             does not point toward the guard set 
%
% Syntax:  
%    [res,R] = checkFlow(obj,guard,R,options)
%
% Inputs:
%    obj - object of class location
%    guard - guard set (class: constrained hyperplane)
%    R - list of reachable sets
%    options - struct containing the algorithm settings
%
% Outputs:
%    res - 1 if flow point toward the guard, 0 otherwise
%    R - updated list of reachable sets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/guardIntersect

% Author:       Niklas Kochdumper
% Written:      23-December-2019             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    sys = obj.contDynamics;

    % adapt the guard set such that the normal vector of the hyperplane points
    % toward the outside of the invariant set
    guard = adaptGuard(obj,guard,R);
    
    % loop over all reachable sets
    R_ = cell(length(R),1);
    counter = 1;
    
    for i = 1:length(R)
            
        % check if the flow points toward the guard set
        res = flowInDirection(sys,guard,R{i},options);
        if res
           R_{counter} = R{i}; 
           counter = counter + 1;
        end    
    end
    
    % check if all sets are empty
    if counter > 1
        res = 1;
        R = R_(1:counter-1);
    else
        res = 0;
        R = [];
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = flowInDirection(sys,guard,R,options)
% check if the flow of the system points in the direction of the guard set
    
    % get hyperplane normal direction
    c = guard.h.c./norm(guard.h.c);

    % project reachable set onto the hyperplane
    Rproj = projectOnHyperplane(guard,R);

    % fast check: check if the flow of the center points towards the guard
    options.u = center(options.U);
    
    fcnHan = getfcn(sys,options);
    flow = fcnHan(0,center(Rproj));
    flow = c'*flow;
    
    if flow > 0
       res = 1;
       return;
    end

    % compute interval enclosure of the reachable set
    B = gramSchmidt(guard.h.c);
    
    if isa(Rproj,'zonoBundle')
       Rproj = Rproj.Z{1}; 
    end
    
    if ~isa(Rproj,'zonotope')
       Rproj = zonotope(Rproj); 
    end
    
    R = B * reduce(B'*Rproj,'pca',1);
    
    % compute flow in hyperplane normal direction
    if isa(sys,'linearSys')
        cen = center(R);
        G = generators(R);
        int = interval(-ones(size(G,2),1),ones(size(G,2),1));
        temp = c'*sys.A;
        flow = temp*cen + (temp*G)*int + ...
                (c'*sys.B)*interval(options.U);
        if ~isempty(sys.c)
           flow = flow + c'*sys.c; 
        end
    else
        tay = taylm(R);
        options.u = taylm(interval(options.U),4,'u');
        fcnHan = getfcn(sys,options);
        flow = fcnHan(0,tay);
        flow = c' * flow;
    end
    
    % check if the flow points in the direction of the guard set
    res = supremum(interval(flow)) > 0;
end

function guard = adaptGuard(obj,guard,R)
% adapt the guard set such that the normal vector of the hyperplane points
% toward the outside of the invariant set

    hs = halfspace(guard.h.c,guard.h.d);
    inv = obj.invariant;
    
    % get center of the first reachable set
    if iscell(R{1})
       c = center(R{1}{1}); 
    else
       c = center(R{1}); 
    end
    
    % check if center is on the correct side of the hyperplane
    if in(inv,c)
       if ~in(hs,c)
           guard = conHyperplane(-guard.h.c,-guard.h.d);
       end
    else
       if in(hs,c)
           guard = conHyperplane(-guard.h.c,-guard.h.d);
       end
    end
end

%------------- END OF CODE --------------