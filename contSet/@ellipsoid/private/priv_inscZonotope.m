function Z = priv_inscZonotope(E,m,mode)
% priv_inscZonotope - inner-approximates a non-degenerate ellipsoid by a
%    zonotope
%
% Syntax:
%    Z = priv_inscZonotope(E,m,mode)
%
% Inputs:
%    E - ellipsoid object
%    m - number of generators
%    mode - computation of zonotope norm (see zonotope/norm_)
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z = zonotope(E,'inner:norm',20);
% 
%    figure; hold on;
%    plot(E); plot(Z);
%
% References:
%    [1] P. Leopardi, "Distributing points on the sphere," School of Mathematics
%        and Statistics. PhD thesis. University of New South Wales, 2007
%    [2] The Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox,
%        https://de.mathworks.com/matlabcentral/fileexchange/13356-eq_sphere_partitions
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_encZonotope

% Authors:       Victor Gassmann, Matthias Althoff, Tobias Ladner
% Written:       18-September-2019
% Last update:   27-January-2021 (MA, degenerate case implemented)
%                05-February-2025 (TL, case n==m and symmetry breaking)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out dimension
n = dim(E);
% extract center
c = center(E);
% compute transformation matrix s.t. T*E == unit hypersphere
Tinv = sqrtm(E.Q);
% compute "uniform" distribution of m points on unit hypersphere, see [1],[2]
if n == 1
    % return exact result
    G = 1;
elseif n == m
    G = eye(n);
elseif mod(m,2) == 0
    % such cases result in aligned generators
    % -> choose twice as many and discard half of it
    G = eq_point_set(n-1,m*2);
    G = G(:,1:m);
else
    G = eq_point_set(n-1,m);
    % check if aligned
    if size(generators(compact(zonotope(zeros(n,1),G),'all')),2) ~= m
        % break symmetry; choose one more and discard it
        % TL: this can probably be improved by enlarging the other
        % generators by the discarded one
        Z = priv_inscZonotope(E,m+1,mode);
        G = Z.G(:,2:end);
    end
end

% init zonotope
Z = zonotope(zeros(n,1),G);

% compute zonotope norm
R = norm_(Z,2,mode);
if isnan(R) 
    R = norm_(Z,2,'ub');
end

% we want the zonotope to be enclosed in the ellipsoid, so we scale
% zonotope(.,G) such that is barely contained in unit hypersphere,
% and apply inverse transform
Z = c + 1/R*Tinv*Z;

% ------------------------------ END OF CODE ------------------------------
