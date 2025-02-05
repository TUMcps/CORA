function [res,cert,scaling] = contains_(C,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if a capsule contains a set or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(C,S,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    C - capsule object
%    S - contSet object or single point
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of C
%       will be detected as lying in C, which can be useful to counteract
%       errors originating from floating point errors.
%    maxEval - Currently has no effect
%    certToggle - if set to 'true', cert will be computed (see below),
%       otherwise cert will be set to NaN.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below), otherwise scaling will be set to inf.
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in C, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in C).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(C - C.center) + C.center contains S.
%
% Example: 
%    C1 = capsule([0;0],[-2;2],2);
%    C2 = capsule([-0.5;0],[-0.5;1],1);
%    C3 = capsule([-2;0],[-0.5;1],1);
%
%    contains(C1,C2)
%    contains(C1,C3)
%
%    figure; hold on;
%    plot(C1,[1,2],'b');
%    plot(C2,[1,2],'g');
%
%    figure; hold on;
%    plot(C1,[1,2],'b');
%    plot(C3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, zonotope/contains_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       20-November-2019
% Last update:   15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------
    
    % init result
    res = true;
    cert = NaN;
    scaling = Inf;
    
    if scalingToggle
        throw(CORAerror('CORA:notSupported',...
                    "The computation of the scaling factor for " + ...
                    "capsules is not yet implemented."));
    end

    % point in capsule containment
    if isnumeric(S)
        
        res = false(1,size(S,2));
        cert = true(1,size(S,2)); % Whatever the result, it is guaranteed
                                  % to be correct here
        
        for i = 1:size(S,2)
            res(i) = aux_containsPoint(C,S(:,i));
        end

        return

    % capsule is empty

    elseif representsa(C, 'emptySet')
        res = representsa(S, 'emptySet');
        cert = true;
        if res
            scaling = 0;
        else
            scaling = inf;
        end

    % empty set is trivially contained
    elseif isa(S,'emptySet') || representsa(S,'emptySet')
        res = true;
        cert = true;
        scaling = 0;
        return

    % fullspace is trivially not contained
    elseif isa(S,'fullspace')
        res = false;
        cert = true;
        scaling = inf;
        return
        
    % capsule in capsule containment    
    elseif isa(S,'capsule')
        
        res = aux_containsCapsule(C,S);
        cert = true;
        return

    % polytopic set in capsule containment    
    elseif isa(S,'conZonotope') || isa(S,'interval') || isa(S,'polytope') ...
            || isa(S,'zonoBundle') || isa(S,'zonotope') || isa(S,'polygon')
        
        % compute vertices
        V = vertices(S);
        
        % check containment for each vertex
        % this evaluation is always exact
        for i = 1:size(V,2)
            if ~aux_containsPoint(C,V(:,i))
                res = false;
                cert = true;
                return;
            end
        end
        cert = true;

        return

        % non polytopic set in capsule containment
    elseif isa(S,'ellipsoid') || isa(S,'taylm') || isa(S,'polyZonotope') ...
            || isa(S,'conPolyZono') || isa(S,'spectraShadow')

        % if the user wants the exact result, throw an error
        if ~strcmp(method, 'approx')
            throw(CORAerror('CORA:noExactAlg',C,S));
        end
       
        % check containment with over-approximating zonotope
        res = contains_(C,zonotope(S),'exact',tol,maxEval,certToggle,scalingToggle);
        % this evaluation is not necessarily exact
        if res
            cert = true;
        else
            cert = false;
        end

        return

    else
        throw(CORAerror('CORA:noops',C,S));
    end
    
end


% Auxiliary functions -----------------------------------------------------

function res = aux_containsCapsule(C1,C2)
% checks if the capsule obj2 is contained in the capsule obj2

    if isempty(C2.g)
       res = aux_containsSphere(C1,C2.c,C2.r); 
    else
       res1 = aux_containsSphere(C1,C2.c+C2.g,C2.r);
       res2 = aux_containsSphere(C1,C2.c-C2.g,C2.r);
       
       res = res1 & res2;
    end
end

function res = aux_containsSphere(C,c,r)
% checks if the capsule contains the hypersphere defined by the center 
% and the radius
    
    % check case where capsule is a hypersphere (no generator)
    if isempty(C.g)
        
        tmp = norm(C.c-c) + r;
        res = tmp < C.r | withinTol(tmp,C.r);
       
    else
       
        ng = norm(C.g);
        g_ = C.g/ng;
        proj = (c-C.c)'*g_;
        
        % check if center is in hypercylinder
        if abs(proj) < norm(C.g)
           
            % compute distance to axis
            diff = (C.c + proj*g_) - c;
            
            % check if distance to axis is smaller than the radius
            tmp = norm(diff) + r;
            res = tmp < C.r | withinTol(tmp,C.r);
            
        else    
            % check if point is in upper or lower hypersphere
            tmp = norm(C.c + C.g - c) + r;
            res1 = tmp < C.r | withinTol(tmp,C.r);
            tmp = norm(C.c - C.g - c) + r;
            res2 = tmp < C.r | withinTol(tmp,C.r);
            
            res = res1 | res2;
        end  
    end
end

function res = aux_containsPoint(C,p)
% checks if a point is contained in the capsule
   
    % get object properties
    g = C.g;
    r = C.r;
    c = C.c;

    % check case where capsule is a hyperball
    if isempty(g)
        
       res = aux_inSphere(c,r,p);
       
    else
       
        ng = norm(g);
        g_ = g/ng;
        proj = (p-c)'*g_;
        
        % check if point is in hypercylinder
        if abs(proj) < norm(g)
           
            % compute distance to axis, check if smaller than the radius
            dist = norm((c + proj*g_) - p);
            res = dist < r | withinTol(dist,r);
            
        else    
            % check if point is in upper or lower hypersphere
            res = aux_inSphere(c+g,r,p) | aux_inSphere(c-g,r,p);
        end  
    end
end


function res = aux_inSphere(c,r,p)
% checks if a point is contained in the hypersphere defined by center and
% radius
    
    tmp = norm(p-c);
    res = tmp < r | withinTol(tmp,r);

end

% ------------------------------ END OF CODE ------------------------------
