function [A,B,E] = platoon(n)
% platoon - scalable platooning system from [1, (9)]
%    the state vector is a sequence of concatenated vectors
%        [e_i, dot(e)_i, a_i]  with  i=1...n
%    for n following vehicles, where e_i is the distance between the i-1th
%    and ith vehicle and a_i is the acceleration of the ith vehicle;
%    the control input vector is composed of the individual accelerations
%
% Syntax:  
%    [A,B,E] = platoon(n)
%
% Inputs:
%    n - state dimension
%
% Outputs:
%    A - state matrix
%    B - input matrix
%    E - disturbance matrix
%
% Reference:
%    [1] I. Ben Makhlouf, S. Kowalewski. "Optimizing Safe Control of a
%        Networked Platoon of Trucks Using Reachability", ARCH14-15,
%        pp. 169-179.

% Author:       Mark Wetzlinger
% Written:      07-August-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% system represented as x' = Ax + Bu
% however, CORA currently does not support E matrices

% gamma = 1/T_i, where T_i the time constant of the drivetrain of vehicle i
gamma = 0.5;

% state matrix:
diagpart = kron(eye(n),[0 1 0; 0 0 -1; 0 0 -gamma]);
offdiagpart = blkdiag(diag(repmat([0;0;1],n-1,1),-2),0);
A = diagpart + offdiagpart;

% input matrix: acceleration for each vehicle
B = kron(eye(n),[0;0;gamma]);

% disturbance matrix: acceleration of leading vehicle
E = zeros(3*n,1);
E(2) = 1;

%------------- END OF CODE --------------
