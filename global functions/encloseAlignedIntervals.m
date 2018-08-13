function Z = encloseAlignedIntervals(Zall,D,guard)
% encloseAlignedIntervals - encloses a set of aligned intervals with a
%                           zonotope bundle
%
% Syntax:  
%    Z = encloseAlignedIntervals(Zall,D)
%
% Inputs:
%    Zall - cell array containing the aligned intervals
%    D - cell array containing the corresponding orthogonal basis
%    guard - guard set that is active
%
% Outputs:
%    Z - enclosing zonotope bundle
%
% Example: 
%
% Other m-files required: not specified
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    Znew = cell(length(D),1);

    % loop over all different basis
    for k = 1:length(D)

        infi = inf * ones(size(D{k},1),1);
        sup = -inf * ones(size(D{k},1),1);

        % loop over all reachable sets that intersect the guard
        for i = 1:length(Zall)
           if ~iscell(Zall{i})
               if ~isempty(Zall{i}.intervals{k})
                   infi = min(infi,infimum(Zall{i}.intervals{k}));
                   sup = max(sup,supremum(Zall{i}.intervals{k}));
               end
           else
               for j = 1:length(Zall{i})
                   if ~isempty(Zall{i}{j}.intervals{k})
                       infi = min(infi,infimum(Zall{i}{j}.intervals{k}));
                       sup = max(sup,supremum(Zall{i}{j}.intervals{k}));
                   end
               end
           end
        end
        
        % intersect the resulting interval with the feasible region form
        % the constrained hyperplane guard set
        if isa(guard,'constrainedHyperplane')
            [infi,sup] = intersectInterval(infi,sup,guard,D{k});
        end

        % transform resulting interval back to original space
        temp = interval(infi,sup);       
        Znew{k} = zonotope(D{k}*[mid(temp),diag(rad(temp))]);
    end
    
    % generate resulting zonotopeBundle object
    if length(Znew) > 1
        Z = zonotopeBundle(Znew);
    else
        Z = zonotope(Znew{1}); 
    end
end


% Auxiliary functions -----------------------------------------------------

function [infi,sup] = intersectInterval(infi,sup,cHp,D)
% Remove the parts of the interval that are located outside the feasible
% region defined by all constraints

    % constraint matrices
    C = cHp.C;
    d = cHp.d;
    
    % transform to current space
    C = C*D;
    
    inter = interval(infi,sup);
    m = mid(inter);
    
    % loop over all constraints
    for i = 1:size(C,1)
       
        % check if the constraint hyperplane intersects the interval
        if isIntersecting(halfspace(C(i,:),d(i)),inter)
            
            % Remove interval parts outside the feasible region
            dom = reduceDomConstr(interval(infi,sup),C(i,:),d(i));
            infi = infimum(dom);
            sup = supremum(dom);
            
        elseif C(i,:)*m > d(i)   
            
            % interval not in feasible region -> break
            infi = [];
            sup = [];
            break;
        end       
    end
end


function dom = reduceDomConstr(dom,C,d)
% Remove the parts of the interval that are located outside the feasible
% region defined by one constraint


    domNew = dom;

    % loop over all dimensions
    for i = 1:length(dom)
        if abs(C(i)) > 1e-10
           temp = C;
           temp(i) = 0;
           domNew(i) = (-temp*dom + d)/C(i);
        end
    end

    domNew = domNew & dom;  
    decVec = sign(C);
    
    % reduce the upper or lower bound of the domain 
    for i = 1:length(dom)
        if decVec(i) == 1
           dom(i) = interval(infimum(dom(i)),supremum(domNew(i))); 
        elseif decVec(i) == -1
           dom(i) = interval(infimum(domNew(i)),supremum(dom(i)));
        end
    end   
end

%------------- END OF CODE --------------