function [res,cert,scaling] = priv_zonotopeContainment_ellipsoidSampling(E,Z,tol,N,scalingToggle)
% priv_zonotopeContainment_ellipsoidSampling - Solves the
% ellipsoid-in-zonotope containment problem by using the Shenmaier vertex
% sampling algorithm described in [1].
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_ellipsoidSampling(E,Z,tol,N,scalingToggle)
%
% Inputs:
%    E - ellipsoid object, inbody
%    Z - zonotope object, circumbody
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
%           verified. For example, if res=false and cert=true, E is
%           guaranteed to not be contained in Z, whereas if res=false and
%           cert=false, nothing can be deduced (E could still be
%           contained in Z).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(Z - center(Z)) + center(Z) contains E.
%           For priv_zonotopeContainment_ellipsoidSampling, this is a lower
%           bound.
%           Note that computing this scaling factor may significantly
%           increase the runtime.
%
% Example:
%    Z = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    E = ellipsoid([1 0; 0 2]);
% 
%    % The function priv_zonotopeContainment_ellipsoidSampling is called
%    % implicitly by contains
%    contains(Z,E,'sampling')
% 
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(E,[1,2],'g');
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

Q = E.Q;
q = E.q;

% Computing the projection-representation of E
[U,Sigma,~] = svd(Q);
s = diag(Sigma);
s_plus = 1./s;
s_plus = sqrt(s_plus);
s_plus(abs(s) < 1000 * eps) = 0;

G = U' * diag(s_plus);
c = q;

% Generate points uniformly on the boundary of E
p = randn([dim(Z) N]);
p = p./sqrt(sum(p.^2));

p = G * p + c;

% should scaling be computed?
if scalingToggle
    [res,~,scaling] = contains_(Z,p,'exact',tol,N,true,true);
    
    res = all(res);
    cert = ~res;
    scaling = max(scaling);
else
    % no scaling
    scaling = NaN;
    res = true;
    for i=1:N
        [res,~,~] = contains_(Z,p(:,i),'exact',tol,N,false,false);
        if ~res
            break
        end
    end
    cert = ~res;
end

end

% ------------------------------ END OF CODE ------------------------------
