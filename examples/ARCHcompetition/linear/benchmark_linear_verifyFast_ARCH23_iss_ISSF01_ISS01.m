function text = benchmark_linear_verifyFast_ARCH23_iss_ISSF01_ISS01
% benchmark_linear_verifyFast_ARCH23_iss_ISSF01_ISS01 - iss benchmark from
%     the 2023 ARCH competition
%
% Syntax:
%    text = benchmark_linear_verifyFast_ARCH23_iss_ISSF01_ISS01()
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

% Parameters --------------------------------------------------------------

% initial set
R0 = interval(-0.0001*ones(270,1),0.0001*ones(270,1));
params.R0 = zonotope(R0);

% uncertain inputs
U = interval([0;0.8;0.9],[0.1;1;1]);
params.U = zonotope(U);

% final time
params.tFinal = 20;

options = struct();
options.verifyAlg = 'reachavoid:supportFunc';

% Specification -----------------------------------------------------------

% forall t: -7e-4 <= y3 <= 7e-4
d = 7e-4;
P1 = polytope([0, 0, 1], -d); 
P2 = polytope([0, 0, -1], -d);
spec = specification({P1,P2},'unsafeSet');


% System Dynamics ---------------------------------------------------------

% load system matrices
load iss.mat A B C

% construct the linear system object
sys = linearSys('iss',A,B,[],C);


% Verification ------------------------------------------------------------

% min steps needed: 1240
[res,fals,savedata] = verify(sys,params,options,spec);

disp("specifications verified: " + res);
disp("computation time: " + savedata.tComp);

% Return value ------------------------------------------------------------

text = ['Spacestation,ISSF01-ISS01,',num2str(res),',',num2str(savedata.tComp)];

% ------------------------------ END OF CODE ------------------------------
