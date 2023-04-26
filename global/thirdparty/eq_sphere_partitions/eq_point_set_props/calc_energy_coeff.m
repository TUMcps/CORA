function coeff = calc_energy_coeff(dim,N,s,energy)
%CALC_ENERGY_COEFF Coefficient of second term in expansion of energy
%
%Syntax
% coeff = calc_energy_coeff(d,N,s,energy);
%
%Description
% COEFF = CALC_ENERGY_COEFF(dim,N,s,ENERGY) sets COEFF to be the coefficient of
% the second term of an expansion of ENERGY with the same form as the expansion 
% of E(dim,N,s), the minimum r^(-s) energy of a set of N points on S^dim.
%
% Specifically, for s not equal to 0, COEFF is the solution to
%
% ENERGY == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim),
%
% and for s == 0 (the logarithmic potential), COEFF is the solution to
%
% ENERGY == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The argument ENERGY must an array of real numbers of the same array size as N.
% The result COEFF will be an array of the same size as N.
%
%Notes
% 1) The energy expansion is not valid for N == 1, and in particular,
%
% EQ_ENERGY_COEFF(dim,N,0,energy) := 0.
%
% 2) For s > 0, [KuiS98 (1.6) p524] has
%
% E(dim,N,s) == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim) + ...
% 
% where SPHERE_INT_ENERGY(dim,s) is the energy integral of the r^(-s) potential
% on S^dim.
%
% The case s == 0 (logarithmic potential) can be split into subcases.
% For s == 0 and dim == 1, E(1,N,0) is obtained by equally spaced points on S^1,
% and the formula for the log potential for N equally spaced points on a circle
% gives
%
% E(1,N,0) == (-1/2) N LOG(N) exactly.
%
% For s == 0 and dim == 2, [SafK97 (4) p7] has
%
% E(2,N,0) == (SPHERE_INT_ENERGY(2,0)/2) N^2 + COEFF N LOG(N) + o(N LOG(N)).
%
% In general, for s == 0,
%
% E(dim,N,0) == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N) + ... 
%
% with sphere_int_energy(1,0) == 0.
%
% CALC_ENERGY_COEFF just uses this general formula for s == 0, so for s == 0 and
% dim == 1, the coefficient returned is actually the coefficient of the first
% non-zero term.
%
%Examples
% > N=2:6
%  N =
%       2     3     4     5     6
%  
% > energy=eq_energy_dist(2,N,0)
%  energy =
%     -0.6931   -1.3863   -2.7726   -4.4205   -6.2383
%  
% > calc_energy_coeff(2,N,0,energy)
%  ans =
%     -0.2213   -0.1569   -0.2213   -0.2493   -0.2569
%
%See also
% EQ_ENERGY_DIST, EQ_ENERGY_COEFF

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Compute the energy coefficient: subtract the first term in the expansion of
% the minimum energy and divide by the power of N in the second term.
%
if s > 0
    %
    % The first term in the expansion of the minimum energy.
    % Ref: [KuiS98 (1.6) p524]
    %
    first_term = (sphere_int_energy(dim,s)/2) * N.^2;
    coeff = (energy-first_term) ./ (N.^(1+s/dim));
else
    %
    % Flatten N into a row vector.
    %
    shape = size(N);
    n_partitions = prod(shape);
    N = reshape(N,1,n_partitions);
    %
    % Refs for s==0, dim == 2: 
    % [SafK97 (4) p. 7] [Zho95 (5.6) p. 68, 3.11 - corrected) p. 42]
    %
    first_term = (sphere_int_energy(dim,s)/2) * N.^2;
    %
    % Avoid division by zero.
    %
    coeff = zeros(size(N));
    neq1 = (N ~= 1);
    coeff(neq1) = (energy(neq1)-first_term(neq1)) ./ (N(neq1) .* log(N(neq1)));
    %
    % Reshape output to same array size as original N.
    %
    coeff = reshape(coeff,shape);
    %
end
%
% end function

function energy = sphere_int_energy(dim,s)
%SPHERE_INT_ENERGY Energy integral of r^(-s) potential
%
%Syntax
% energy = sphere_int_energy(d,s);
%
%Description
% ENERGY = SPHERE_INT_ENERGY(dim,s) sets ENERGY to be the energy integral 
% on S^dim of the r^(-s) potential, defined using normalized Lebesgue measure.
%
% Ref for s > 0: [KuiS98 (1.6) p524]
% Ref for s == 0 and dim == 2: SafK97 (4) p. 7] 
% For s == 0 and dim >= 2, integral was obtained using Maple:
%     energy = (1/2)*(omega(dim)/omega(dim+1)* ...
%          int(-log(2*sin(theta/2)*(sin(theta))^(dim-1),theta=0..Pi),
%     where omega(dim+1) == area_of_sphere(dim).

if s ~= 0
    energy = (gamma((dim+1)/2)*gamma(dim-s)/(gamma((dim-s+1)/2)*gamma(dim-s/2)));
elseif dim ~= 1
    energy = (psi(dim)-psi(dim/2)-log(4))/2;
else
    energy = 0;    
end
%
% end function
