function P = normalizeConstraints(P,varargin)
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
%    P - mptPolytope object
%    type - (optional) 'b' (default): normalize offset vectors b and be
%                      'A': normalize norm of constaints in A and Ae to 1
%
% Outputs:
%    P - polytope object
%
% Example: 
%    A = [2 1; -2 3; -2 0; -1 -4; 2 -3; 4 -1];
%    b = [2 3 -1 2 0 1]';
%    P = mptPolytope(A,b);
%    P = normalizeConstraints(P,'A');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      01-December-2022
% Last update:  13-December-2022 (MW, add A normalization)
% Last revision:---

%------------- BEGIN CODE --------------

% set default type
type = setDefaultValues({'b'},varargin);

% check input arguments
inputArgsCheck({{P,'att','mptPolytope'}, ...
                {type,'str',{'b','A'}}});

% read out properties
A = P.P.A; b = P.P.b;
Ae = P.P.Ae; be = P.P.be;

if strcmp(type,'b')
    % normalize offset vectors b and be to -1|0|1

    % normalize inequality constraints
    if ~isempty(A)
        % inequality constraints where A(:,i)*x <= 0 are left unchanged
    
        % find indices of constraints with b > 0 and b < 0
        idx_plus = b > 0;
        idx_neg = b < 0;
    
        % divide constraints with b > 0 by b
        A(idx_plus,:) = A(idx_plus,:) ./ b(idx_plus);
        b(idx_plus) = 1;
    
        % divide constraints with b < 0 by |b|
        A(idx_neg,:) = A(idx_neg,:) ./ abs(b(idx_neg));
        b(idx_neg) = -1;
    
    end
    
    % normalize equality constraints
    if ~isempty(Ae)
        % equality constraints where Ae(:,i)*x = 0 are left unchanged
    
        % find indices of constraints with be > 0 and be < 0
        idx_plus = be > 0;
        idx_neg = be < 0;
    
        % divide constraints with be > 0 by be
        Ae(idx_plus,:) = Ae(idx_plus,:) ./ be(idx_plus);
        be(idx_plus) = 1;
    
        % divide constraints with be < 0 by be
        Ae(idx_neg,:) = Ae(idx_neg,:) ./ be(idx_neg);
        be(idx_neg) = 1;
    
    end

elseif strcmp(type,'A')
    % normalize norm of constraints in A and Ae to 1

    % normalize inequality constraints
    if ~isempty(A)
        normA = vecnorm(A',2,1);
        A = (A' ./ normA)';
        b = b ./ normA';
    end
    
    % normalize equality constraints
    if ~isempty(Ae)    
        normA = vecnorm(Ae',2,1);
        Ae = (Ae' ./ normA)';
        be = be ./ normA';
    end

end

% init resulting polytope
S = struct('A',A,'b',b,'Ae',Ae,'be',be);
P = mptPolytope(S);

%------------- END OF CODE --------------