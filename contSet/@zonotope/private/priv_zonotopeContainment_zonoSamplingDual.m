function [res,cert,scaling] = priv_zonotopeContainment_zonoSamplingDual(Z1,Z2,tol,N,scalingToggle)
% priv_zonotopeContainment_zonoSamplingDual - Solves the zonotope containment problem by using
%    the Shenmaier halfspace sampling algorithm described in [1].
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_zonoSamplingDual(Z1,Z2,tol,N,scalingToggle)
%
% Inputs:
%    Z1 - zonotope object, inbody
%    Z2 - zonotope object, circumbody
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of Z
%       will be detected as lying in Z, which can be useful to counteract
%       errors originating from floating point errors.
%    N - Number of random samples
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
%           For priv_zonotopeContainment_zonoSampling, this is a lower bound.
%           Note that computing this scaling factor may significantly
%           increase the runtime.
%
% Example:
%    Z1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    Z2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    Z3 = Z2 + [3;0];
% 
%    % The function priv_zonotopeContainment_zonoSamplingDual is called implicitly by contains
%    contains(Z1,Z2,'samplingDual')
%    contains(Z1,Z3,'samplingDual')
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
%    
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
%
% References:
%    [1] Kulmburg A., Brkan I., Althoff M.,: Search-based and Stochastic
%        Solutions to the Zonotope and Ellipsotope Containment Problems
%        (to appear)
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

% init values
c = center(Z1);
G = Z1.generators;

d = center(Z2);
H = Z2.generators;
n = size(H,1);
ell = size(H, 2);

G_prime = [G c-d];

% compute scaling
scaling = 0;
for i = 1:N
    % Compute a random facet of the zonotope
    permutation = randperm(ell);
    ind = permutation(1:n-1);
    Q = H(:,ind);
    lambda = ndimCross(Q);
    
    b = sum(abs(lambda' * H));
    
    % If b==0, the origin lies on the facet, indicating that the zonotope
    % is degenerate. This case can be dealt with manually, by interpreting
    % lambda as 0
    if b == 0
        lambda = 0.*lambda;
    else
        lambda = lambda / b;
    end
    
    s = sum(abs(G_prime' * lambda));
    
    scaling = max([scaling s]);
    
    if ~scalingToggle && scaling > 1 + tol
        break
    end
end

% check if scaling is within tol
res = scaling <= 1+tol;
cert = ~res;

end

% ------------------------------ END OF CODE ------------------------------
