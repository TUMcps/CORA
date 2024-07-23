function P_out = normalizeConstraints(P,varargin)
% normalizeConstraints - normalizes all inequality/equality constraints;
%    type 'b' (default) normalizes the entries in the offset vector b, so
%    that we have
%       A(i,:) x <= 1
%       A(i,:) x <= 0
%       A(i,:) x <= -1
%    for the inequality constraints and
%       Ae(i,:) x = 1
%       Ae(i,:) x = 0
%    for the equality constraints; alternatively, using type 'A', the norm
%    of the constraints in A and Ae is normalized to 1.
%
% Syntax:
%    P = normalizeConstraints(P)
%    P = normalizeConstraints(P,type)
%
% Inputs:
%    P - polytope object
%    type - (optional) 'b', 'be' (default): normalize offset vectors b and be
%                      'A', 'Ae': normalize norm of constaints in A and Ae to 1
%
% Outputs:
%    P_out - polytope object
%
% Example: 
%    A = [2 1; -2 3; -2 0; -1 -4; 2 -3; 4 -1];
%    b = [2 3 -1 2 0 1]';
%    P = polytope(A,b);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       01-December-2022
% Last update:   13-December-2022 (MW, add A normalization)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% can only normalize constraints if constraints are given
if ~P.isHRep.val
    throw(CORAerror('CORA:specialError',...
        'Halfspace representation not computed... nothing to normalize!'));
end

% set default type
type = setDefaultValues({'b'},varargin);

% check input arguments
inputArgsCheck({{P,'att','polytope'}, ...
                {type,'str',{'b','be','A','Ae'}}});

% init output polytope
P_out = polytope(P);

if strcmp(type,'b') || strcmp(type,'be')
    % normalize offset vectors b and be to -1|0|1

    % normalize inequality constraints
    if ~isempty(P.A_.val)
        % inequality constraints where A(:,i)*x <= 0 are left unchanged
    
        % find indices of constraints with b > 0 and b < 0
        idx_plus = P.b_.val > 0;
        idx_neg = P.b_.val < 0;
    
        % divide constraints with b > 0 by b
        if any(idx_plus)
            P_out.A_.val(idx_plus,:) = P.A_.val(idx_plus,:) ./ P.b_.val(idx_plus);
            P_out.b_.val(idx_plus) = 1;
        end
    
        % divide constraints with b < 0 by |b|
        if any(idx_neg)
            P_out.A_.val(idx_neg,:) = P.A_.val(idx_neg,:) ./ abs(P.b_.val(idx_neg));
            P_out.b_.val(idx_neg) = -1;
        end
    
    end
    
    % normalize equality constraints
    if ~isempty(P.Ae_.val)
        % equality constraints where Ae(:,i)*x = 0 are left unchanged
    
        % find indices of constraints with be > 0 and be < 0
        idx_nonzero = P.be_.val > 0 | P.be_.val < 0;
    
        % divide constraints with be =!= 0 by be
        if any(idx_nonzero)
            P_out.Ae_.val(idx_nonzero,:) = P.Ae_.val(idx_nonzero,:) ./ P.be_.val(idx_nonzero);
            P_out.be_.val(idx_nonzero) = 1;
        end
        
    end

elseif strcmp(type,'A') || strcmp(type,'Ae')
    % normalize norm of constraints in A and Ae to 1
    % skip constraints of the form 0*x <= ... or 0*x == ...

    % normalize inequality constraints
    if ~isempty(P.A_.val)
        normA = vecnorm(P.A_.val',2,1);
        idx_nonzero = ~withinTol(normA,0);
        if any(idx_nonzero)
            P_out.A_.val(idx_nonzero,:) = (P.A_.val(idx_nonzero,:)' ./ normA(idx_nonzero))';
            P_out.b_.val(idx_nonzero) = P.b_.val(idx_nonzero) ./ normA(idx_nonzero)';
        end
    end
    
    % normalize equality constraints
    if ~isempty(P.Ae_.val)    
        normA = vecnorm(P.Ae_.val',2,1);
        idx_nonzero = ~withinTol(normA,0);
        if any(idx_nonzero)
            P_out.Ae_.val(idx_nonzero,:) = (P.Ae_.val(idx_nonzero,:)' ./ normA(idx_nonzero))';
            P_out.be_.val(idx_nonzero) = P.be_.val(idx_nonzero) ./ normA(idx_nonzero)';
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
