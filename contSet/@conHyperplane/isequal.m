function res = isequal(hyp1,hyp2,varargin)
% isequal - checks if two constrained hyperplanes are equal
%
% Syntax:
%    res = isequal(hyp1,hyp2)
%    res = isequal(hyp1,hyp2,tol)
%
% Inputs:
%    hyp1 - conHyperplane object
%    hyp2 - conHyperplane object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp1 = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    hyp2 = hyp1;
%    hyp3 = conHyperplane(halfspace([1;-1],0),[1 0;-1 0],[2;2]);
%
%    isequal(hyp1,hyp2)
%    isequal(hyp1,hyp3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-September-2019
% Last update:   06-June-2022
%                01-March-2024 (MW, update deprecated syntax)
% Last revision: 10-January-2023 (MW, expand overly restrictive comparison)

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% set default value
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{hyp1,'att','conHyperplane'};
                {hyp2,'att',{'conHyperplane','polytope'}}; ...
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% if second object is a polytope, check conversion to conHyperplane
[res,hyp2] = representsa_(hyp2,'conHyperplane',tol);
if ~res; return; end

% assume true
res = true;

% 1. check if normal vector of halfspace is (anti-)parallel

% normalization: a*x = b -> a'*b = 1 (or 0)
[a_norm_1,b_norm_1] = aux_normalizeEquality(hyp1.a,hyp1.b);
[a_norm_2,b_norm_2] = aux_normalizeEquality(hyp2.a,hyp2.b);

if ~all(withinTol([a_norm_1;b_norm_1],[a_norm_2;b_norm_2],tol))
    res = false; return
end

% 2. check if constraints are equal

% ensure that constraints are given
hyp1_constrained = ~isempty(hyp1.C);
hyp2_constrained = ~isempty(hyp2.C);
if hyp1_constrained && hyp2_constrained
    % normalize inequalities
    [a_norm_1,b_norm_1] = aux_normalizeInequalities(hyp1.C,hyp1.d);
    [a_norm_2,b_norm_2] = aux_normalizeInequalities(hyp2.C,hyp2.d);

    % compare inequality constraints
    if ~compareMatrices([a_norm_1 b_norm_1]',[a_norm_2 b_norm_2]',tol)
        res = false; return
    end
elseif xor(hyp1_constrained,hyp2_constrained)
    % caution: there may be special case, where the constraints have no
    % effect, then the hyperplanes would still be equal...
    res = false; return
end

end


% Auxiliary functions -----------------------------------------------------

function [a_norm,b_norm] = aux_normalizeEquality(a,b)
% normalize a halfspace a*x = b to either a*x = 1 or a*x = 0

if withinTol(b,0)
    % a*x = 0 ... normalize a to unit length
    a_norm = a' / vecnorm(a);
    % fix first non-zero entry to be positive
    idx_nonzero = find(~withinTol(a_norm,0),1,'first');
    if a_norm(idx_nonzero) < 0
        a_norm = -1*a_norm;
    end
    b_norm = 0;
else
    % a*x = b ... normalize to a*x = 1
    a_norm = a' / b;
    b_norm = 1;
end

end

function [C_norm,d_norm] = aux_normalizeInequalities(C,d)
% normalize a halfspace C*x <= d such that
%   either C(i,:)*x <= 1
%          C(i,:)*x <= -1,
%       or C(i,:)*x <= 0 (with C(i,:) having unit length

% init output variables
C_norm = zeros(size(C));
d_norm = zeros(size(d));

% loop over each inequality
for i=1:length(d_norm)

    if withinTol(d(i),0)
        % C(i,:)*x <= 0 ... normalize C(i,:) to unit length
        C_norm(i,:) = C(i,:) / vecnorm(C(i,:));
        % fix first non-zero entry to be positive
        idx_nonzero = find(~withinTol(C_norm(i,:),0),1,'first');
        if C_norm(i,idx_nonzero) < 0
            C_norm(i,:) = -1*C_norm(i,:);
        end
        d_norm(i) = 0;

    else
        % C(i,:)*x <= d ... normalize to C(i,:)*x = -1/+1
        C_norm(i,:) = C(i,:) / abs(d(i));
        d_norm(i) = sign(d(i)) * 1;
    end

end

end

% ------------------------------ END OF CODE ------------------------------
