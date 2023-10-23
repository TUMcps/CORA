function res = test_linearSysDT_reach_01_5dim()
% test_linearSysDT_reach_01_5dim - unit test function of linear reachability 
%    discrete-time analysis with uncertain inputs
%
% Checks the solution of the linearSys class for a 5-dimensional example;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = test_linearSysDT_reach_01_5dim()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff
% Written:       24-March-2020
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

% system matrix
A = [0.95 -0.15 0 0 0;...
    0.15 0.95 0 0 0; ...
    0 0 0.9 0.05 0; ...
    0 0 -0.05 0.9 0; ...
    0 0 0 0 0.92];
dim_x = length(A);

% input matrix
B = 1;

% sampling time
dt = 0.04;

fiveDimSys = linearSysDT('fiveDimSys',A,B,dt); %initialize system


% Parameters --------------------------------------------------------------

params.tFinal = 5; %final time
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x)); %initial set
params.U = zonotope(zeros(dim_x,1),0.02*diag([0.1, 0.3, 0.1, 0.3, 0.3])); %input


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder=200; %zonotope order


% Reachability Analysis (zonotope) ----------------------------------------

R = reach(fiveDimSys, params, options);

% enclose reachable set by interval
IH = interval(R.timePoint.set{end-1});

%saved result
IH_saved = interval( ...
           [-0.1283222186855075; -0.1231839187600670; -0.0409825827512041; -0.0585969262587082; -0.0749684720882265], ...
           [0.1331532444946742; 0.1451671122132107; 0.0409897108679144; 0.0585982732675874; 0.0750331447277620]);

%final result
res_zono = isequal(IH,IH_saved,1e-8);


% Reachability Analysis (zonotope bundles) --------------------------------

%generate zonotope bundle
Z0{1} = params.R0;
Z0{2} = params.R0 + [-0.1; 0; 0.1; 0; 0];
params.R0 = zonoBundle(Z0);


R = reach(fiveDimSys, params, options);

% enclose reachable set by interval
IH = interval(R.timePoint.set{end-1});

%saved result
IH_saved = interval( ...
    [-0.1283222186855075; -0.1231839187600670; -0.0409823708730643; -0.0585969262587082; -0.0749684720882265], ...
    [0.1324828890131164; 0.1447383080221113; 0.0409897108679144; 0.0585981287398916; 0.0750331447277620]);

%final result
res_zonoBundles = isequal(IH,IH_saved,1e-8);

%result of different set representations
res = res_zono && res_zonoBundles;

% ------------------------------ END OF CODE ------------------------------
