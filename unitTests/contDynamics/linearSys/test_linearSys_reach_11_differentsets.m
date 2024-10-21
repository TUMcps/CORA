function res = test_linearSys_reach_11_differentsets
% test_linearSys_reach_11_differentsets - unit test for reach,
%    where many possible different sets (U,W,V) are checked for completion;
%
% Syntax:
%    res = test_linearSys_reach_11_differentsets
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
% Written:       19-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define linearSys --------------------------------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
states = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
inputs = size(B,2);

% constant offset: n x 1
c = 0.05 * [-4; 2; 3; 1];

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];
outputs = 2;

% feedthrough matrix: q x m
D = [0 0 1;
     0 0 0];

% constant input: q x 1
k = [0; 0.02];

% disturbance matrix: n x r
E = [1 0.5; 0 -0.5; 1 -1; 0 1];
dists = size(E,2);

% noise matrix: q x s
F = [1; 0.5];
noises = size(F,2);


% initialize different linearSys-objects
sys_A = linearSys(A,1);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);
sys_ABcCDkE = linearSys(A,B,c,C,D,k,E);
sys_ABcCDkEF = linearSys(A,B,c,C,D,k,E,F);
list = {sys_A, sys_AB, sys_ABC, sys_ABCD, sys_ABcCDk, sys_ABcCDkE, sys_ABcCDkEF};


% reachability settings and model parameters ------------------------------

options.timeStep = 0.05;
options.taylorTerms = 6;
options.zonotopeOrder = 50;

% time horizon
params.tFinal = 0.20;
dt_steps = params.tFinal / options.timeStep;
% initial set
params.R0 = zonotope(10*ones(states,1),0.5*diag(ones(states,1)));
% input vector
u_n = randn(states,dt_steps+1);
u_n_noD = randn(states,dt_steps);
u_m = randn(inputs,dt_steps+1);
u_m_noD = randn(inputs,dt_steps);
% input set
U_n = interval([-4;-2;-1;-1],[-2;2;1;0]);
U_m = interval([-2;-1;-1],[2;1;0]);
% disturbance set
W = zonotope(0.02+zeros(dists,1),0.02*diag(ones(dists,1)));
% sensor noise set
V = zonotope(-0.01+zeros(noises,1),0.01*diag(ones(noises,1)));


% reach -------------------------------------------------------------------

methods = {'standard','fromStart','wrapping-free'}; % 'adaptive'

for j = 1:length(list)
    
    sys = list{j};
    
    for k = 1:length(methods)
        
        if k > 1
            options_cont.linAlg = methods{k};
        end
    
        % no input set, disturbance, or sensor noise
        R = reach(sys,params,options);

        % input set
        if sys.nrOfInputs == states
            params.U = U_n;
        elseif sys.nrOfInputs == inputs
            params.U = U_m;
        end
        R = reach(sys,params,options);

        % input vector
        if sys.nrOfInputs == states
            if ~isempty(sys.D)
                params.u = u_n;
            else
                params.u = u_n_noD;
            end
        elseif sys.nrOfInputs == inputs
            if ~isempty(sys.D)
                params.u = u_m;
            else
                params.u = u_m_noD;
            end
        end
        R = reach(sys,params,options);

        % disturbance
        if sys.nrOfDisturbances == dists
            params.W = W;
        end
        R = reach(sys,params,options);

        % sensor noise
        if sys.nrOfNoises == noises
            params.V = V;
        end
        R = reach(sys,params,options);

        % remove fields for clean next iteration
        params = rmiffield(params,'U');
        params = rmiffield(params,'u');
        params = rmiffield(params,'W');
        params = rmiffield(params,'V');
        
        if k > 1
            options_cont = rmfield(options_cont,'linAlg');
        end
    end
        
end

% all checks ok
res = true;

end

% ------------------------------ END OF CODE ------------------------------
