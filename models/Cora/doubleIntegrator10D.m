function [A,B,E] = doubleIntegrator10D
% doubleIntegrator10D - double integrator system with uncontrollable
%    subspace from [1, (18)]
%
% Syntax:
%    [A,B,E] = doubleIntegrator10D
%
% Inputs:
%    -
%
% Outputs:
%    A - state matrix
%    B - input matrix
%    E - disturbance matrix
%
% Reference:
%    [1] L. Yang and N. Ozay, "Scalable Zonotopic Under-Approximation of
%        Backward Reachable Sets for Uncertain Linear Systems", IEEE
%        Control Systems Letters (6), 2022.

% Authors:       Mark Wetzlinger
% Written:       08-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system represented as x' = Ax + Bu
% however, CORA currently does not support E matrices

% state matrix
A = [0 1 0 0 0 0  1     0  0  1;
     0 0 0 0 0 0  0     0  0  0;
     0 0 0 1 0 0  0    -1  0  0;
     0 0 0 0 0 0  0     0  0  0;
     0 0 0 0 0 1  0     0  1  0;
     0 0 0 0 0 0  0     0  0  0;
     0 0 0 0 0 0 -1e-2  1  0  0;
     0 0 0 0 0 0 -1e-2 -1  0  0;
     0 0 0 0 0 0 -1e-4  0  0  2;
     0 0 0 0 0 0  0     0 -2 -1e-4];

% input matrix
B = zeros(10,3);
B(2,1) = 1;
B(4,2) = 1;
B(6,3) = 1;

% disturbance matrix
E = eye(10);

% ------------------------------ END OF CODE ------------------------------
