function sys = linearSysDT(sys,dt)
% linearSysDT - convert linear continuous-time system to discrete-time
%               system
%
% Syntax:
%    sys = linearSysDT(sys,dt)
%
% Description:
%    Converts a continuous-time linear system to a equivalent discrete-time 
%    linear system according to Eq. (6) in [1]
%
% Inputs:
%    sys - continuous-time linear system (class: linearSys)
%    dt - sampling time
%
% Outputs:
%    sys - discrete-time linear system (class: linearSysDT)
%
% Example:
%    % continuous-time system
%    A = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
%    B = [0 0;0 0;1 0;0 1];
%    C = [1 0 0 0;0 1 0 0];
%
%    sys = linearSys(A,B,[],C);
% 
%    % convert to discrete-time system
%    Ts = 0.1;
%    sys_ = linearSysDT(sys,Ts);
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert system matrix A_ = e^A*dt
A = expm(sys.A*dt);

% convert input matrix B_ = A^-1 * (e^A*dt - I) * B
temp1 = eye(size(sys.A,1)) * dt;
temp2 = temp1;
cnt = 2;

while true
    temp1 = temp1 * dt/cnt * sys.A;
    temp2 = temp2 + temp1;
    cnt = cnt + 1;
    if all(all(abs(temp1) < eps)) || cnt > 1000
        break;
    end
end

B = temp2 * sys.B;

% convert constant input c_ = A^-1 * (e^A*dt - I) * c
c = zeros(sys.dim,1);
if ~isempty(sys.c)
    c = temp2 * sys.c;
end

% construct resulting discrete time system
sys = linearSysDT(A,B,c,sys.C,sys.D,sys.k,dt);
    
end

% ------------------------------ END OF CODE ------------------------------
