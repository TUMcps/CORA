function [res,ind] = check(obj,S)
% check - checks if a set satisfies the specification
%
% Syntax:  
%    [res,ind] = check(obj,R)
%
% Inputs:
%    obj - specification object
%    S - contSet object
%
% Outputs:
%    res - 1 if set satisfies the specification, 0 if not
%    ind - index of the specification that is violated
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Author:       Niklas Kochdumper
% Written:      29-May-2020             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    ind = [];

    % loop over all specifications
    for i = 1:size(obj,1)
       
        % different types of specifications
        switch obj(i,1).type
             
            case 'invariant'
                res = checkInvariant(obj(i,1).set,S);
                
            case 'unsafeSet'
                res = checkUnsafeSet(obj(i,1).set,S);
                
            case 'safeSet'
                res = checkSafeSet(obj(i,1).set,S);
                
            case 'custom'
                res = checkCustom(obj(i,1).set,S);
        end
        
        % return as soon as one specification is violated
        if ~res
            ind = i;
            return;
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = checkUnsafeSet(set,S)
% check if reachable set intersects the unsafe sets

    if iscell(S)
        res = 1;
        for i = 1:length(S)
           res = ~isIntersecting(set,S{i}); 
           if ~res
              return; 
           end
        end   
    else
        res = ~isIntersecting(set,S);
    end
end

function res = checkSafeSet(set,S)
% check if reachable set is inside the safe set

    if iscell(S)
        res = 1;
        for i = 1:length(S)
           res = in(set,S{i}); 
           if ~res
              return; 
           end
        end   
    else
        res = in(set,S);
    end
end

function res = checkCustom(func,S)
% check if the reachable set satisfies a user provided specification

    if iscell(S)
        res = 0;
        for i = 1:length(S)
            res = func(S{i});
            if res
               return; 
            end
        end
    else
        res = func(S);
    end
end

function res = checkInvariant(set,S)
% check if reachable set intersects the unsafe sets

    if iscell(S)
        res = 0;
        for i = 1:length(S)
            res = isIntersecting(set,S{i},'approx');
            if res
               return; 
            end
        end
    else
        res = isIntersecting(set,S,'approx');
    end
end

%------------- END OF CODE --------------