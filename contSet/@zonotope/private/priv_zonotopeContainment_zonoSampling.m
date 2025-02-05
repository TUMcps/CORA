function [res,cert,scaling] = priv_zonotopeContainment_zonoSampling(Z1,Z2,tol,N,scalingToggle)
% priv_zonotopeContainment_zonoSampling - Solves the zonotope containment problem by using the
%    Shenmaier vertex sampling algorithm described in [1].
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_zonoSampling(Z1,Z2,tol,N,scalingToggle)
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
%    % The function priv_zonotopeContainment_zonoSampling is called implicitly by contains
%    contains(Z1,Z2,'sampling')
%    contains(Z1,Z3,'sampling')
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

G = generators(Z1);
c = center(Z1);

alphas = sign(rand(size(G,2),N)-0.5);

% should scaling be computed?
if scalingToggle
    [res,~,scaling] = contains_(Z2,G*alphas+c,'exact',tol,N,true,true);
    
    res = all(res);
    cert = ~res;
    scaling = max(scaling);
    return
else
    % no scaling
    scaling = NaN;
    res = true;
    for i=1:N
        p = G * alphas(:,i) + c;
        [res,~,~] = contains_(Z2,p,'exact',tol,N,false,false);
        if ~res
            break
        end
    end
    cert = ~res;
    return
end

end

% ------------------------------ END OF CODE ------------------------------
