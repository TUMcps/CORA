function [A,B,E] = quadrotor
% quadrotor - quadrotor system from [1], linearized about the hover
%    condition
%
% Syntax:  
%    [A,B,E] = platoon(n)
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
%    [1] S. Kaynama, C.J. Tomlin. "Benchmark: Flight Envelope Protection in
%        Autonomous Quadrotors", ARCH14-15.
%    [2] F. Gruber and M. Althoff. "Scalable Robust Safety Filter with
%        Unknown Disturbance Bounds", TAC.

% Author:       Mark Wetzlinger
% Written:      07-August-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% system represented as x' = Ax + Bu
% however, CORA currently does not support E matrices

% constant
g = 9.81;

% state matrix [1, Appendix A]:
% x1, x2, x3: spatial position
% x4, x5, x5: spatial velocity
% x7, x8, x9: angular positions
% x10, x11, x12: angular velocities
A = zeros(12);
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,8) = -g;
A(5,7) = g;
A(7,10) = 1;
A(8,11) = 1;
A(9,12) = 1;

% input matrix [1, Appendix A]:
% u1: total normalized thrust
% u2, u3, u4: second-order derivatives of angular positions
B = zeros(12,4);
B(6,1) = 1;
B(10,2) = 1;
B(11,3) = 1;
B(12,4) = 1;

% disturbance matrix: [2, Sec. V-D]:
% w1, w2, w3: affects spatial velocities
E = zeros(12,3);
E(4,1) = 1;
E(5,2) = 1;
E(6,3) = 1;

%------------- END OF CODE --------------

