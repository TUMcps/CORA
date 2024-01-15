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
    if ~isempty(P.A)
        % inequality constraints where A(:,i)*x <= 0 are left unchanged
    
        % find indices of constraints with b > 0 and b < 0
        idx_plus = P.b > 0;
        idx_neg = P.b < 0;
    
        % divide constraints with b > 0 by b
        P_out.A(idx_plus,:) = P.A(idx_plus,:) ./ P.b(idx_plus);
        P_out.b(idx_plus) = 1;
    
        % divide constraints with b < 0 by |b|
        P_out.A(idx_neg,:) = P.A(idx_neg,:) ./ abs(P.b(idx_neg));
        P_out.b(idx_neg) = -1;
    
    end
    
    % normalize equality constraints
    if ~isempty(P.Ae)
        % equality constraints where Ae(:,i)*x = 0 are left unchanged
    
        % find indices of constraints with be > 0 and be < 0
        idx_nonzero = P.be > 0 | P.be < 0;
    
        % divide constraints with be =!= 0 by be
        if any(idx_nonzero)
            P_out.Ae(idx_nonzero,:) = P.Ae(idx_nonzero,:) ./ P.be(idx_nonzero);
            P_out.be(idx_nonzero) = 1;
        end
        
    end

elseif strcmp(type,'A') || strcmp(type,'Ae')
    % normalize norm of constraints in A and Ae to 1
    % skip constraints of the form 0*x <= ... or 0*x == ...

    % normalize inequality constraints
    if ~isempty(P.A)
        normA = vecnorm(P.A',2,1);
        idx_nonzero = ~withinTol(normA,0);
        if any(idx_nonzero)
            P_out.A(idx_nonzero,:) = (P.A(idx_nonzero,:)' ./ normA(idx_nonzero))';
            P_out.b(idx_nonzero) = P.b(idx_nonzero) ./ normA(idx_nonzero)';
        end
    end
    
    % normalize equality constraints
    if ~isempty(P.Ae)    
        normA = vecnorm(P.Ae',2,1);
        idx_nonzero = ~withinTol(normA,0);
        if any(idx_nonzero)
            P_out.Ae(idx_nonzero,:) = (P.Ae(idx_nonzero,:)' ./ normA(idx_nonzero))';
            P_out.be(idx_nonzero) = P.be(idx_nonzero) ./ normA(idx_nonzero)';
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
