function text = benchmark_linear_verifyFast_ARCH23_beam_CBC01
% benchmark_linear_verifyFast_ARCH23_beam_CBC01 - beam benchmark from the
%     2023 ARCH competition
%
% Syntax:
%    text = benchmark_linear_verifyFast_ARCH23_beam_CBC01()
%
% Inputs:
%    -
%
% Outputs:
%    text - char array

% Authors:       Mark Wetzlinger
% Written:       23-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Model Derivation --------------------------------------------------------

% nodes in model
N = 100;
% node of interest
node = round(0.7*N);

% constants
rho = 7.3e-4;   % density
L = 200;        % length of beam
Q = 1;          % cross-section area (renamed from A)
E = 30e6;       % Young's modulus

ell = L/N;      % length of individual discrete element

% mass matrix (NxN)
M = (rho*Q*ell) / 2 * diag([2*ones(N-1,1);1]);
Minv = M^(-1);

% load
F = zonotope(10000,100);

% tridiagonal matrix (NxN)
mat = zeros(N);
mat(1,1) = 2; mat(1,2) = -1;
mat(N,N-1) = -1; mat(N,N) = 1;
for r=2:N-1
    mat(r,1+(r-2)) = -1;
    mat(r,2+(r-2)) = 2;
    mat(r,3+(r-2)) = -1;
end
% stiffness matrix (NxN)
K = E*Q/ell * mat;
% damping matrix (NxN)
a = 1e-6;
b = 1e-6;
D = a*K + b*M;

% state matrix (damped)
A = [zeros(N) eye(N); -Minv*K -Minv*D];

% nr of states
dim_x = length(A);

% constant inputs
A_C = [ A, [zeros(dim_x-1,1);1]; zeros(1,dim_x+1) ];
B_C = zeros(size(A_C,1),1);


% Parameters --------------------------------------------------------------

params.tFinal = 0.01;

% initial set: bar at rest
params.R0 = cartProd( zonotope(zeros(dim_x,1)), Minv(end,end)*F );

% input set
params.U = zonotope(0);

options = struct();
options.verifyAlg = 'reachavoid:supportFunc';

% System Dynamics ---------------------------------------------------------

C = zeros(1,2*N+1);
C(1,2*node) = 1;

% construct linear system objects
sys = linearSys('beam',A_C,B_C,[],C);


% Specification -----------------------------------------------------------

% forall t: y1 <= 74
spec = specification(halfspace(1,74),'safeSet');


% Verification ------------------------------------------------------------

% min steps needed: 695
[res,fals,savedata] = verify(sys,params,options,spec);


% Return value ------------------------------------------------------------

text = ['Beam,CBC01,',num2str(res),',',num2str(savedata.tComp)];

% ------------------------------ END OF CODE ------------------------------
