function res = test_linearSys_simulateRandom_01
% test_linearSys_simulateRandom_01 - unit test for simulation of linearSys,
%    where many possible different setting are checked for completion;
%    note: the numerical correctness of the result is not (yet) checked!
%
% Syntax:
%    res = test_linearSys_simulateRandom_01
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

tFinal = 2;
params = {};
% case: only R0
params{1}.tFinal = tFinal;
params{1}.R0 = zonotope([ones(5,1),0.1*diag(ones(5,1))]);
% case: R0 and U
params{2}.tFinal = tFinal;
params{2}.R0 = params{1}.R0;
params{2}.U = zonotope(interval([0.9; -0.25; -0.1], ...
                             [1.1; 0.25; 0.1]));


% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
dim_x = length(A);          % n
dim_u = dim(params{2}.U);   % m
dim_y = 2;                  % p
% case: A, B \in \R{n x m}
B = randn(dim_x,dim_u);
sys{1} = linearSys('fiveDimSys',A,B);
% case: A, B, c \in \R{n}
c = randn(dim_x,1); 
sys{2} = linearSys('fiveDimSys',A,B,c);
% case: A, B, c, C \in \R{p x n}
C = randn(dim_y,dim_x);
sys{3} = linearSys('fiveDimSys',A,B,c,C);
% case: A, B, c, C, D \in \R{p x m}
D = randn(dim_y,dim_u);
sys{4} = linearSys('fiveDimSys',A,B,c,C,D);
% case: A, B, c, C, D, k \in \R{p}
k = rand(dim_y,1);
sys{5} = linearSys('fiveDimSys',A,B,c,C,D,k);


% Random Simulations ------------------------------------------------------

simOpt = {};
% case: only one point
simOpt{1}.points = 1;
simOpt{1}.fracVert = 0.5;
simOpt{1}.fracInpVert = 0.5;
simOpt{1}.nrConstInp = 10;
% case: no input changes
simOpt{2}.points = 1;
simOpt{2}.fracVert = 0.5;
simOpt{2}.fracInpVert = 0.5;
simOpt{2}.nrConstInp = 1;
% case: only extreme points
simOpt{3}.points = 1;
simOpt{3}.fracVert = 1;
simOpt{3}.fracInpVert = 1;
simOpt{3}.nrConstInp = 10;


% Check for completion ----------------------------------------------------

res = [];
% loop over all defined systems
for s=1:length(sys)
    % loop over all defined model parameters
    for p=1:length(params)
        % loop over all settings
        for o=1:length(simOpt)
            try
                simRes = simulateRandom(sys{s}, params{p}, simOpt{o});
                res(end+1,1) = true;
            catch
                % run-time error
                res(end+1,1) = false;
            end
        end
    end
end

% combine results
res = all(res);


end

% ------------------------------ END OF CODE ------------------------------
