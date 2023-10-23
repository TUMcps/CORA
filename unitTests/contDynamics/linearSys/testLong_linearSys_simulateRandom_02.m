function res = testLong_linearSys_simulateRandom_02
% testLong_linearSys_simulateRandom_02 - unit test for simulation of linearSys,
%    where all possible set representations for the initial set are used
%    note: the numerical correctness of the result is not (yet) checked!
%
% Syntax:
%    res = testLong_linearSys_simulateRandom_02
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;


% Model Parameters --------------------------------------------------------

params.tFinal = 5;
% R0 using all set representations
R0{1} = interval(0.9*ones(5,1),1.1*ones(5,1));          % interval
R0{2,1} = capsule(R0{1});                               % capsule
R0{3} = conPolyZono(R0{1});                             % conPolyZono
R0{4} = conZonotope(R0{1});                             % conZonotope
R0{5} = ellipsoid(R0{1});                               % ellipsoid
R0{6} = polytope(R0{1});                                % polytope
R0{7} = polyZonotope(R0{1});                            % polyZonotope
R0{8} = zonotope(R0{1});                                % zonotope
R0{9} = zonoBundle({R0{8},R0{8}+0.1*ones(5,1)});        % zonoBundle
% U using all set representations
U{1} = interval([0.9; -0.25; -0.1],[1.1; 0.25; 0.1]);   % interval
U{2,1} = capsule(U{1});                                 % capsule
U{3} = conPolyZono(U{1});                               % conPolyZono
U{4} = conZonotope(U{1});                               % conZonotope
U{5} = ellipsoid(U{1});                                 % ellipsoid
U{6} = polytope(U{1});                               % polytope
U{7} = polyZonotope(U{1});                              % polyZonotope
U{8} = zonotope(U{1});                                  % zonotope
U{9} = zonoBundle({U{8},U{8}+0.1*ones(3,1)});           % zonoBundle


% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
dim_x = length(A);
dim_u = dim(U{1});
B = randn(dim_x,dim_u);
c = randn(dim_x,1);
% output equation
dim_y = 2;
C = randn(dim_y,dim_x);
D = randn(dim_y,dim_u);
k = rand(dim_y,1);
% most general case: A, B, c, C, D, k
sys = linearSys('fiveDimSys',A,B,c,C,D,k);


% Check for completion ----------------------------------------------------

% loop over all R0 / U
res = false(length(R0),1);
for r=1:length(R0)
    params.R0 = R0{r};
    params.U = U{r};
    try
        simRes = simulateRandom(sys, params);
        res(r,1) = true;
    catch ME
        % run-time error
        res(r,1) = false;
    end
end

% combine results
res = all(res);


end

% ------------------------------ END OF CODE ------------------------------
