function coeffs = findBernsteinPoly(f,l,u,n)
% findBernsteinPoly - finds a polynomial approximating f on [l,u] using
%     bernstein polynomials
%
% Syntax:
%    coeffs = findBernsteinPoly(f,l,u,n)
%
% Inputs:
%    f - function handle to approximate
%    l,u - bounds of domain
%    n - order of the polynomial
%
% Outputs:
%    coeffs - coefficients of the polynomial
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       31-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% scale domain to [0,1] and compute polynomial according to
% https://en.wikipedia.org/wiki/Bernstein_polynomial#Approximating_continuous_functions

% scale domain
f_norm = @(x) f(x*(u-l) + l);

% compute bernstein polynomial
coeffs_norm = zeros(1,n+1);
for v=0:n
    coeffs_norm = coeffs_norm + f_norm(v/n) * aux_b(v,n);
end

% transform polynomial back to [l,u] domain
coeffs = zeros(1,n+1);
P = 1; % Pascal's triangle
for i=0:n
    coeffs(end-i:end) = coeffs(end-i:end) + ...
        coeffs_norm(end-i) .* P .* (l/(u-l)).^(0:i) .* (1/(u-l)).^(i:-1:0);
    
    % prepare for next iteration
    P = [1 P(2:end)+P(1:end-1) 1];
end

end


% Auxiliary functions -----------------------------------------------------

function coeffs = aux_b(v,n)
    % returns the bernstein polynomial b_{v,n}(x) in normal form
    coeffs = zeros(1,n+1);
    for l=v:n
        coeffs(end-l) = nchoosek(n,l)*nchoosek(l,v)*(-1)^(l-v);
    end

end

% ------------------------------ END OF CODE ------------------------------
