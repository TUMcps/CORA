function I = interval(hs)
% interval - converts a halfspace to an interval object (exact conversion
%    only possible if normal vector of halfspace is axis-aligned)
%
% Syntax:
%    I = interval(hs)
%
% Inputs:
%    hs - halfspace object
%
% Outputs:
%    I - interval object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
if representsa_(hs,'emptySet',eps)
    I = interval.empty(dim(hs)); return
end

% dimension
n = dim(hs);

% check which dimension of normal vector is (close to) zero
zerodims = withinTol(hs.c,0,eps);

if nnz(zerodims) >= n-1
    % exact conversion possible: normal vector is axis-aligned
    
    if all(zerodims)
        % unbounded -> should be prevented by constructor
        throw(CORAerror('CORA:specialError'));

    else
        % init lower and upper bounds
        lb = -Inf(n,1);
        ub = Inf(n,1);

        % check sign of value at non-zero dimension
        idx = find(~zerodims,1,'first');
        if hs.c(idx) < 0
            % update lower bound
            lb(idx) = hs.d / hs.c(idx);

        else
            % update upper bound
            ub(idx) = hs.d / hs.c(idx);

        end

        % instantiate interval
        I = interval(lb,ub);

    end    

else
    % convert to polytope and use that function
    I = interval(polytope(hs));

end

% ------------------------------ END OF CODE ------------------------------
