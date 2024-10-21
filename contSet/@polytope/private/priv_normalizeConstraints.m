function [A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,type)
% priv_normalizeConstraints - normalizes the constraints
%
% Syntax:
%    [A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,type)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    type - (optional) 'b', 'be': normalize offset vectors b and be
%                      'A', 'Ae': normalize norm of constraints in A and Ae to 1
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(type,'b') || strcmp(type,'be')
    % normalize offset vectors b and be to -1|0|1

    % normalize inequality constraints
    if ~isempty(A)
        % inequality constraints where A(:,i)*x <= 0 are left unchanged
    
        % find indices of constraints with b > 0 and b < 0
        idx_plus = b > 0;
        idx_neg = b < 0;
    
        % divide constraints with b > 0 by b
        if any(idx_plus)
            A(idx_plus,:) = A(idx_plus,:) ./ b(idx_plus);
            b(idx_plus) = 1;
        end
    
        % divide constraints with b < 0 by |b|
        if any(idx_neg)
            A(idx_neg,:) = A(idx_neg,:) ./ abs(b(idx_neg));
            b(idx_neg) = -1;
        end
    
    end
    
    % normalize equality constraints
    if ~isempty(Ae)
        % equality constraints where Ae(:,i)*x = 0 are left unchanged
    
        % find indices of constraints with be > 0 and be < 0
        idx_nonzero = be > 0 | be < 0;
    
        % divide constraints with be =!= 0 by be
        if any(idx_nonzero)
            Ae(idx_nonzero,:) = Ae(idx_nonzero,:) ./ be(idx_nonzero);
            be(idx_nonzero) = 1;
        end
        
    end

elseif strcmp(type,'A') || strcmp(type,'Ae')
    % normalize norm of constraints in A and Ae to 1
    % skip constraints of the form 0*x <= ... or 0*x == ...

    % normalize inequality constraints
    if ~isempty(A)
        normA = vecnorm(A',2,1);
        idx_nonzero = ~withinTol(normA,0);
        if any(idx_nonzero)
            A(idx_nonzero,:) = (A(idx_nonzero,:)' ./ normA(idx_nonzero))';
            b(idx_nonzero) = b(idx_nonzero) ./ normA(idx_nonzero)';
        end
    end
    
    % normalize equality constraints
    if ~isempty(Ae)    
        normA = vecnorm(Ae',2,1);
        idx_nonzero = ~withinTol(normA,0);
        if any(idx_nonzero)
            Ae(idx_nonzero,:) = (Ae(idx_nonzero,:)' ./ normA(idx_nonzero))';
            be(idx_nonzero) = be(idx_nonzero) ./ normA(idx_nonzero)';
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
