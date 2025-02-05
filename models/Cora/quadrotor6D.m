function [A,B,E] = quadrotor6D(xlin,ulin)
% quadrotor6D - linearized quadrotor system from [1, (42)]
%
% Syntax:
%    [A,B,E] = quadrotor6D(xlin,ulin)
%
% Inputs:
%    xlin - six-dimensional linearization point (state)
%    ulin - two-dimensional linearization point (input)
%
% Outputs:
%    A - state matrix
%    B - input matrix
%    E - disturbance matrix
%
% Reference:
%    [1] I.M. Mitchell, J. Budzis, A. Bolyachevets. "Invariant, Viability
%        and Discriminating Kernel Under-Approximation via Zonotope
%        Scaling". arXiv, 2019.
%    [2] A. Sasfi, M. Zeilinger, J. Koehler. "Robust adaptive MPC using
%        control contraction metrics", Automatica, 2023.

% Authors:       Mark Wetzlinger
% Written:       17-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system represented as x' = Ax + Bu + Ew
% however, CORA currently does not support E matrices

% % [2], but re-ordered: [1, 2, 4, 5, 3, 6]
% % constants
% L = 0.25;
% m = 0.486;
% J = 0.00383;
% g = 9.81;
% 
% % state matrix [2, Sec. 4]
% A = zeros(6);
% A(1,3) = cos(xlin(5));
% A(1,4) = -sin(xlin(5));
% A(2,3) = sin(xlin(5));
% A(2,4) = cos(xlin(5));
% A(3,5) = 0;
% A(3,6) = 0;
% A(4,5) = 0;
% A(4,6) = 0;
% A(5,6) = 1;
% 
% 
% % input matrix [2, Sec. 4]
% B = zeros(6,2);
% B(4,1) = 1/m;
% B(4,2) = 1/m;
% B(6,1) = l/J;
% B(6,2) = -l/J;
% 
% 
% % disturbance matrix [2, Sec. 4]
% E = zeros(6,1);
% E(4,1) = 0;
% E(5,1) = 0;


% [1]
% constants
K = 0.89/1.4;
g = 9.81;
d0 = 70;
d1 = 17;
n0 = 55;

% state matrix [1, (42)]:
% x1, x2: horizontal/vertical position
% x3, x4: horizontal/vertical velocity
% x5, x6: roll and roll velocity
A = zeros(6);
A(1,3) = 1;
A(2,4) = 1;
A(3,5) = K * ulin(1) * cos(xlin(5));
A(4,5) = -K * ulin(1) * sin(xlin(5));
A(5,6) = 1;
A(6,5) = -d0;
A(6,6) = -d1;

% input matrix [1, (42)]:
% u1: total thrust
% u2: desired roll angle
B = zeros(6,2);
B(3,1) = K * sin(xlin(5));
B(4,1) = K * cos(xlin(5));
B(6,2) = n0;

% disturbance matrix: [1, (42)]:
E = zeros(6,2);
E(3,1) = 1;
E(4,2) = 1;

% drift term
% w(3,1) = -K * ulin(1) * sin(xlin(5));
% w(4,1) = K * ulin(1) * cos(xlin(5)) - g;

% ------------------------------ END OF CODE ------------------------------
