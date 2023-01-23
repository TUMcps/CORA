function [res,hyp] = isConHyperplane(obj)
% isConHyperplane - check if a polytope can be equivalently represented as
%                   a constrained hyperplane
%
% Syntax:  
%    [res,hyp] = isConHyperplane(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    res - true/false whether equivalently representable as a constrained
%          hyperplane
%    hyp - equivalent conHyperplane object
%
% Example: 
%    A = [1 1;1 0;-1 -1];
%    b = [1;2;-1];
%    P = mptPolytope(A,b);
%
%    [res,hyp] = isConHyperplane(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane

% Author:       Niklas Kochdumper
% Written:      18-May-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = false;
    hyp = [];

    % check if equality constraints exist
    if ~isempty(obj.P.Ae)
        
        % first equality constraint defines hyperplane
        c = obj.P.Ae(1,:);
        d = obj.P.be(1);
        
        % convert all other equality constraints to inequality constraints
        A = [obj.P.A; obj.P.Ae(2:end,:); -obj.P.Ae(2:end,:)];
        b = [obj.P.b; obj.P.be(2:end); -obj.P.be(2:end)];
        
        % construct constrained hyperplane
        res = true;
        if isempty(A) && isempty(b)
            hyp = conHyperplane(c,d);
        else
            hyp = conHyperplane(c,d,A,b);
        end
        
        
    else
        
        % remove redundant halfspaces
        obj.P = minHRep(obj.P);
        
        % normalize the inequality constraints
        A = obj.P.A;
        b = obj.P.b;
        
        n = sum(A.^2,2);
        A = diag(1./n)*A;
        b = b./n;
        
        % check for the case A*x <= b && A*x >= b => A*x == b
        for i = 1:size(A,1)-1
            
            % check if normal vectors are identical
            [result,idx] = ismembertol(-A(i,:),A(i+1:end,:),'ByRows',true);
            idx = idx(1) + i;
         
            if result
                
                % check if offsets are identical
                if b(i) == -b(idx)
                   
                    % equality constraint detected -> construct hyperplane
                    c = A(i,:);
                    d = b(i);
                    
                    A([i,idx],:) = [];
                    b([i,idx],:) = [];
                    
                    res = true;
                    if isempty(A) && isempty(b)
                        hyp = conHyperplane(c,d);
                    else
                        hyp = conHyperplane(c,d,A,b);
                    end
                    return;
                end
            end
        end
    end
end

%------------- END OF CODE --------------