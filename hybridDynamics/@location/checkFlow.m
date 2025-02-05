function [res,R] = checkFlow(loc,guard,R,params)
% checkFlow - remove all intersections for which the flow of the system
%    does not point toward the guard set (see Sec. 4.4.4 in [1])
%
% Syntax:
%    [res,R] = checkFlow(loc,guard,R,params)
%
% Inputs:
%    loc - location object
%    guard - guard set (class: polytope, levelSet)
%    R - cell-array of reachable sets
%    params - model parameters
%
% Outputs:
%    res - true/false whether flow points toward the guard for any of the
%          intersecting reachable sets
%    R - updated cell-array of reachable sets
%
% References:
%    [1] N. Kochdumper. "Extensions to Polynomial Zonotopes and their
%        Application to Verifcation of Cyber-Physical Systems", PhD Thesis 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/guardIntersect

% Authors:       Niklas Kochdumper
% Written:       23-December-2019             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % read out continuous dynamics
    sys = loc.contDynamics;
    
    % adapt the guard set such that the normal vector of the hyperplane
    % points toward the outside of the invariant set
    if isa(guard,'polytope') && representsa_(guard,'conHyperplane',1e-12)
        guard = aux_adaptGuard(loc,guard,R); 
    elseif isa(guard,'levelSet')
        outside = aux_getOutside(loc,guard,R);
    else
        throw(CORAerror('CORA:notSupported',...
            'Requires levelSet or polytope representing constrained hyperplane.'));
    end
    
    % loop over all reachable sets
    R_ = cell(length(R),1);
    counter = 1;
    
    for i=1:length(R)
            
        % check if the flow points toward the guard set
        if isa(guard,'polytope') && representsa_(guard,'conHyperplane',1e-12)
            res = aux_flowInDirection(sys,guard,R{i},params);
        elseif isa(guard,'levelSet')
            res = aux_flowInDirectionLevelSet(sys,guard,R{i},outside,params);
        end
        
        if res
           R_{counter} = R{i}; 
           counter = counter + 1;
        end    
    end
    
    % check if all sets are empty
    if counter > 1
        res = true;
        R = R_(1:counter-1);
    else
        res = false;
        R = [];
    end
end


% Auxiliary functions -----------------------------------------------------

% TODO: move to contDynamics?

function res = aux_flowInDirection(sys,guard,R,params)
% check if the flow of the system points in the direction of the guard set
    
    % get hyperplane normal direction
    c = guard.Ae ./ norm(guard.Ae);

    % project reachable set onto the hyperplane
    Rproj = projectOnHyperplane(R,guard);

    % fast check: check if the flow of the center points towards the guard
    params.u = center(params.U);
    params.w = center(params.W);
    
    fcnHan = getfcn(sys,params);
    flow = fcnHan(0,center(Rproj));
    flow = c*flow;
    
    if flow > 0
        res = true;
        return;
    end

    % compute interval enclosure of the reachable set
    B = gramSchmidt(guard.Ae');
    
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
        I = interval(-ones(size(G,2),1),ones(size(G,2),1));
        cA = c*sys.A;
        flow = cA*cen + (cA*G)*I + (c*sys.B)*interval(params.U) + sys.c;
    else
        tay = taylm(R);
        params.u = taylm(interval(params.U),4,'u');
        fcnHan = getfcn(sys,params);
        flow = fcnHan(0,tay);
        flow = c * flow;
    end
    
    % check if the flow points in the direction of the guard set
    res = supremum(interval(flow)) > 0;
end

function res = aux_flowInDirectionLevelSet(sys,guard,R,outside,params)
% check if the flow of the system points in the direction of the guard set

    % fast check: check if the flow of the center points towards the guard
    params.u = center(params.U);
    params.w = center(params.W);
    
    c = center(R);
    fcnHan = getfcn(sys,params);
    flow = outside*guard.der.grad(c)'*fcnHan(0,c);
    
    if flow > 0
       res = true;
       return;
    end

    % compute interval enclosure of the reachable set
    if ~isa(R,'zonotope')
       R = zonotope(R); 
    end
    
    R = reduce(R,'pca',1);
    
    % compute flow in hyperplane normal direction
    tay = taylm(R);
    params.u = taylm(interval(params.U),4,'u');
    fcnHan = getfcn(sys,params);
    flow = outside*(guard.der.grad(tay))'*fcnHan(0,tay);
    
    % check if the flow points in the direction of the guard set
    res = supremum(interval(flow)) > 0;
end

function guard = aux_adaptGuard(loc,guard,R)
% adapt the guard set such that the normal vector of the hyperplane
% points toward the outside of the invariant set

    P = polytope(guard.Ae,guard.be);
    inv = loc.invariant;
    
    % get center of the first reachable set
    if iscell(R{1})
        c = center(R{1}{1}); 
    else
        c = center(R{1}); 
    end
    
    % check if center is on the correct side of the hyperplane
    if contains_(inv,c,'exact',100*eps,0,false,false)
        if ~contains_(P,c,'exact',100*eps,0,false,false)
            guard = polytope([],[],-guard.Ae,-guard.be);
        end
    else
        if contains_(P,c,'exact',100*eps,0,false,false)
            guard = polytope([],[],-guard.Ae,-guard.be);
        end
    end
end

function outside = aux_getOutside(loc,guard,R)
% determine which side of the guard set is the outside of the invariant set
% g(x) >= 0 is outside => outside = 1
% g(x) < 0 is outside => outside = -1

    inv = loc.invariant;
    outside = 1;
    
    % get center of the first reachable set
    if iscell(R{1})
        c = center(R{1}{1}); 
    else
        c = center(R{1}); 
    end
    
    % check if center is on the correct side of the hyperplane
    if contains_(inv,c,'exact',100*eps,0,false,false)
        if guard.funHan(c) > 0
            outside = -1;
        end
    else
        if guard.funHan(c) <= 0
            outside = -1;
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
