function [res,cert,scaling] = priv_zonotopeContainment_pointContainment(Z, p, method, tol, scalingToggle)
% priv_zonotopeContainment_pointContainment - Checks whether a point (or point cloud) is
%    contained in a zonotope.
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_pointContainment(Z, p, method, tol, scalingToggle)
%
% Inputs:
%    Z - zonotope object
%    p - point, or matrix of points
%    method - method used for the containment check.
%       The available options are:
%           - 'exact': Checks for containment by using either
%               'exact:venum' or 'exact:polymax', depending on the number 
%               of generators of Z and the size of p.
%           - 'exact:venum': Checks for containment by evaluating the 
%              Z-norm of p (see Algorithm 1 in [1]).
%           - 'exact:polymax': Checks for containment by maximizing the
%               polyhedral norm w.r.t. Z over p (see Algorithm 2 in [1]).
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of Z
%       will be detected as lying in Z, which can be useful to counteract
%       errors originating from floating point errors.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below).
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, p is
%           guaranteed to not be contained in Z, whereas if res=false and
%           cert=false, nothing can be deduced (p could still be
%           contained in Z).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(Z - center(Z)) + center(Z) contains p.
%           For priv_zonotopeContainment_pointContainment, this is an exact value.
%           Note that computing this scaling factor may significantly
%           increase the runtime.
%
% Example:
%    Z = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    p = [1;0];
%
%    % The function priv_zonotopeContainment_pointContainment will implicitly be called by
%    % contains
%    contains(Z,p)
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(p(1),p(2),'rx');
%
%
% References:
%    [1] A. Kulmburg, M. Althoff.: On the co-NP-Completeness of the
%        Zonotope Containment Problem, European Journal of Control 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/contains_

% Authors:       Adrian Kulmburg
% Written:       05-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

switch method
    % There are basically two algorithms to choose from: Either we
    % compute it by first calculating all the halfspaces of the
    % zonotope, or we solve it using linear programming.
    % The first method is similar to what polymax does, while the
    % second method is similar to venum, so if the user chooses
    % specifically one of these algorithms we choose the corresponding
    % method, i.e., if the user chose polymax, then we need to compute
    % the halfspace representation of Z:
    case 'exact:polymax'
        % This is just the same as transforming Z into a polytope, and
        % doing the containment check there
        P = polytope(Z);
        [res,cert,scaling] = contains_(P,p,'exact',tol,0,true,scalingToggle);
    case 'exact:venum'
        res = false(1,size(p,2));
        cert = NaN(1,size(p,2));
        scaling = inf(1,size(p,2));
        % check for each point whether zonotope norm is <= 1
        for i = 1:size(p,2)
            scaling(i) = zonotopeNorm(Z,p(:,i)-Z.c);
        end
        res = scaling <= 1 | withinTol(scaling,1,tol);
        cert = true;
    case {'exact', 'approx', 'approx:st'}
        % If the user has not specifically chosen between venum and
        % polymax, we need to make the choice by estimating the runtime
        % of each algorithm.

        % For polymax, the relevant quantity to look for is the maximal
        % number of facets of the zonotope.
        % If the zonotope is non-degenerate (which we may assume if
        % #generators >= dim), then this is given as follows:
        if size(Z.generators,2) >= size(Z.generators,1)
            numFacets = 2*nchoosek(size(Z.generators, 2), size(Z.generators,1)-1);
        else
            % If the zonotope is degenerate, we need to replace the dim by
            % the rank
            numFacets = 2*nchoosek(size(Z.generators,2),rank(Z)-1);
        end
        % We can now estimate the approximate runtime of the halfspace
        % method:
        runtime_halfspaceMethod = numFacets * size(Z.generators,1)^4;
        % The runtime of the linprog method on the other hand mainly
        % depends on the number of generators, but also on the number of
        % points we are evaluating this on:
        runtime_LP = (size(Z.generators, 2)+1)^(3.5) * size(p,2);

        % We can now compare runtimes. Additionally, we need to
        % double-check that the halfspace method does not lead to an
        % explosion of halfspaces. Currently, the maximal number of
        % halfspaces allowed is hardcoded in such a way, that
        % numberOfFacets < 100000
        % This has to be taken with a grain of salt, and the user may want
        % to change that value to whatever value he/she desires.
        % In practice, this is only relevant for dimensions <= 3.

        % halfspace conversion always preferrable if the zonotope is in
        % fact a parallelotope
        % check if the halfspaces have already been computed:
        useLP = numFacets > 100000 && runtime_halfspaceMethod > runtime_LP;
        if representsa_(Z,'parallelotope',eps) || ~useLP
            % If all the conditions are met, we call the function again
            % using polymax
            [res,cert,scaling] = priv_zonotopeContainment_pointContainment(Z, p, 'exact:polymax', tol, scalingToggle);
        else
            % Otherwise, we divert to venum
            [res,cert,scaling] = priv_zonotopeContainment_pointContainment(Z, p, 'exact:venum', tol, scalingToggle);
        end
    otherwise
        throw(CORAerror('CORA:noSpecificAlg',method,Z,p));
end

end

% ------------------------------ END OF CODE ------------------------------
