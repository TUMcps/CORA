function res = test_Krylov_crossChecks(~)
% test_Krylov_crossChecks - unit test that cross checks certain aspects of
%     the Krylov method. Further information can be found in [1].
%
% Syntax:
%    res = test_Krylov_crossChecks(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Saad, Y. Analysis of Some Krylov Subspace Approximations to the Matrix 
%        Exponential Operator, SIAM Journal on Numerical Analysis, 1992, 29, 209-228

% Authors:       Matthias Althoff
% Written:       15-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% example 1; Saad 1992
N = 100;
for i=1:N
    A(i,i)=(i+1)/(N+1);
end

%obtain dimension
dim = length(A);
redDim = ceil(0.5*dim);

%set options --------------------------------------------------------------
options.timeStep=1; %time step size for reachable set computation
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
linDyn=linearSys('KrylovTest',A,1); %initialize quadratic dynamics
%--------------------------------------------------------------------------

% use Arnoldi iteration
v = ones(dim,1);
[V,H] = arnoldi(A,v,redDim);

% test H
H_diff = H - V'*A*V;
res(1) = all(all(abs(H_diff)<=1e-7)); 

% initial value solution
v_next = expm(A)*v;
v_next_approx = norm(v)*V*expm(H)*[1;zeros(redDim-1,1)];
error = max(abs(v_next - v_next_approx));
res(2) = (error <= 1e-14);

% alternative
v_next_approx = V*expm(H)*V'*v;
error = max(abs(v_next - v_next_approx));
res(3) = (error <= 1e-6);

% test if Krylov method is a projection approach
I_diff = eye(redDim) - V'*V;
res(4) = all(all(abs(I_diff)<=1e-7));

% all tests passed?
res = all(res);

% ------------------------------ END OF CODE ------------------------------
