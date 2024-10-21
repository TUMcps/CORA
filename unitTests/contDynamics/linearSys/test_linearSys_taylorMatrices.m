function res = test_linearSys_taylorMatrices
% test_linearSys_taylorMatrices - unit test for the computation of
%    auxiliary interval matrices in the computation of reachable sets for
%    linear systems
%
% Syntax:
%    res = test_linearSys_taylorMatrices
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init dynamics (only state matrix relevant)
A = [-1 -4; 4 -1];
sys = linearSys(A);

% set time step size and truncation order
timeStep = 0.1;
truncationOrder = 6;
[E,F,G] = taylorMatrices(sys,timeStep,truncationOrder);

% Taylor series expansion of e^Adt = sum_i (A*dt)^i/i! until truncation
% order plus computed remainder term should contain result of built-in expm
% function (assuming that result is exact)
eAdt_approx = eye(2);
for i=1:truncationOrder
    eAdt_approx = eAdt_approx + (A*timeStep)^i/factorial(i);
end
eAdt_enclosure = eAdt_approx + E;
eAdt_exact = expm(A*timeStep);
% convert to 2*nx1 vector for containment check
assert(contains(eAdt_enclosure(:),eAdt_exact(:)));

% correction matrix for input should be approximately equal inverse of A
% times the correction matrix for the state (due to non-zero remainder term
% we use containment, instead)
invA_F = inv(A)*F;
assert(contains(invA_F(:),G(:)));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
