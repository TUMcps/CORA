function [res,cert,scaling] = priv_zonotopeContainment_DIRECTMaximization(Z1, Z2, tol, maxEval, scalingToggle)
% priv_zonotopeContainment_DIRECTMaximization - Solves the zonotope containment problem by
%    checking whether the maximum value of the Z2-norm at one of the
%    vertices of Z1 exceeds 1+tol, using DIRECT optimization, in case
%    ga optimization is not available.
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_DIRECTMaximization(Z1, Z2, tol, maxEval, scalingToggle)
%
% Inputs:
%    Z1 - zonotope object, inbody
%    Z2 - zonotope object, circumbody
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of Z
%       will be detected as lying in Z, which can be useful to counteract
%       errors originating from floating point errors.
%    maxEval - Number of maximal function evaluations.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below).
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, Z1 is
%           guaranteed to not be contained in Z2, whereas if res=false and
%           cert=false, nothing can be deduced (Z1 could still be
%           contained in Z2).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(Z2 - center(Z2)) + center(Z2) contains Z1.
%           For priv_zonotopeContainment_DIRECTMaximization, this is a lower bound.
%           Note that computing this scaling factor may significantly 
%           increase the runtime.
%
% Example: 
%    Z1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    Z2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    Z3 = Z2 + [3;0];
%    
%    % The function priv_zonotopeContainment_DIRECTMaximization will implicitly be called by
%    % contains
%    contains(Z1,Z2,'opt')
%    contains(Z1,Z3,'opt')
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
%    
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
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

% Retrieve the generator matrix of Z1 and determine its size
G = Z1.G;
m = size(G, 2);

% Prepare objective function to be minimized (see [1]).
Problem.f = @(nu) -zonotopeNorm(Z2, G*nu + Z1.c-Z2.c);

bounds = [-ones(m,1) ones(m,1)];
opts.maxevals = maxEval;
opts.showits = 0;
if scalingToggle
    opts.limit = -inf;
else
    opts.limit = -1-tol;
end
opts.maxits    = inf;
                                        
% Launch the optimization. Note that the fourth argument ensures that the
% points that are tested are integer-points, since a maximum can only
% happen on one of the vertices of Z1, meaning that nu should only have
% the values +-1, i.e., integer values.

[scaling,~,~] = Direct(Problem,bounds,opts);

scaling = abs(scaling);

if scaling > 1 + tol
    res = false;
    cert = true;
else
    res = true;
    cert = false;
end
    
end

% ------------------------------ END OF CODE ------------------------------
