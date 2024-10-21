function linsysDT = linearSysDT(linsys,dt)
% linearSysDT - convert a linear continuous-time system to equivalent
%    discrete-time linear system according to Eq. (6) in [1]
%
% Syntax:
%    linsysDT = linearSysDT(linsys,dt)
%
% Inputs:
%    linsys - linearSys object
%    dt - sampling time
%
% Outputs:
%    linsysDT - linearSysDT object
%
% Example:
%    % continuous-time system
%    A = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
%    B = [0 0;0 0;1 0;0 1];
%    C = [1 0 0 0;0 1 0 0];
%    linsys = linearSys(A,B,[],C);
% 
%    % convert to discrete-time system
%    dt = 0.1;
%    linsysDT = linearSysDT(linsys,dt);
%
% References:
%    [1] B. Schuermann et. al: Guaranteeing Constraints of Disturbed
%        Nonlinear Systems Using Set-BasedOptimal Control in Generator 
%        Space, IFAC 2017
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys

% Authors:       Niklas Kochdumper
% Written:       21-November-2020 
% Last update:   19-November-2021 (MW, minor fixes)
%                02-September-2024 (MW, extend to disturbance/noise matrices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert system matrix A_ = e^A*dt
A_dt = expm(linsys.A*dt);

% convert input matrix
%    B_dt = A^-1 * (e^A*dt - I) * B
%    B_dt = T * B
term_j = eye(linsys.nrOfStates) * dt;
T = term_j;

% maximum a thousand terms
for j=2:1000
    term_j = term_j * dt/j * linsys.A;
    T = T + term_j;
    if all(all(abs(term_j) < eps))
        break;
    end
end

B_dt = T * linsys.B;

% same conversion process for disturbance matrix
E_dt = T * linsys.E;

% convert constant input c_ = A^-1 * (e^A*dt - I) * c
c_dt = T * linsys.c;

% construct resulting discrete time system
linsysDT = linearSysDT(A_dt,B_dt,c_dt,...
                    linsys.C,linsys.D,linsys.k,...
                    E_dt,linsys.F,dt);

% ------------------------------ END OF CODE ------------------------------
