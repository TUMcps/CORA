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

% Author:       Mark Wetzlinger
% Written:      19-November-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

% define linearSys --------------------------------------------------------

% stable system matrix: n x n
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
n = length(A);

% input matrix: n x m
B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
m = size(B,2);

% constant offset: n x 1
c = 0.05 * [-4; 2; 3; 1];

% output matrix: q x n
C = [1 1 0 0;
     0 -0.5 0.5 0];
y = 2;

% throughput matrix: q x m
D = [0 0 1;
     0 0 0];

% constant input: q x 1
k = [0; 0.02];

% initialize different linearSys-objects
sys_A = linearSys(A,1);
sys_AB = linearSys(A,B);
sys_ABC = linearSys(A,B,[],C);
sys_ABCD = linearSys(A,B,[],C,D);
sys_ABcCDk = linearSys(A,B,c,C,D,k);
list = {sys_A, sys_AB, sys_ABC, sys_ABCD, sys_ABcCDk};


% reachability settings and model parameters ------------------------------

options.timeStep = 0.05;
options.taylorTerms = 6;
options.zonotopeOrder = 50;

% time horizon
params.tFinal = 0.20;
dt_steps = params.tFinal / options.timeStep;
% initial set
params.R0 = zonotope(10*ones(n,1),0.5*diag(ones(n,1)));
% input vector
u_n = randn(n,dt_steps+1);
u_n_noD = randn(n,dt_steps);
u_m = randn(m,dt_steps+1);
u_m_noD = randn(m,dt_steps);
% input set
U_n = interval([-4;-2;-1;-1],[-2;2;1;0]);
U_m = interval([-2;-1;-1],[2;1;0]);
% disturbance set
W = zonotope(0.02+zeros(n,1),0.02*diag(ones(n,1)));
% sensor noise set
V_n = zonotope(zeros(n,1),diag(ones(n,1)));
V_y = zonotope(-0.01+zeros(y,1),0.01*diag(ones(y,1)));


% reach -------------------------------------------------------------------

methods = {'standard','fromStart','wrapping-free'}; % 'adaptive'

for j = 1:length(list)
    
    sys = list{j};
    
    for k = 1:length(methods)
        
        if k > 1
            options_cont.linAlg = methods{k};
        end
    
        % no input set, disturbace, or sensor noise
        R = reach(sys,params,options);

        % input set
        if sys.nrOfInputs == n
            params.U = U_n;
        elseif sys.nrOfInputs == m
            params.U = U_m;
        end
        R = reach(sys,params,options);

        % input vector
        if sys.nrOfInputs == n
            if ~isempty(sys.D)
                params.u = u_n;
            else
                params.u = u_n_noD;
            end
        elseif sys.nrOfInputs == m
            if ~isempty(sys.D)
                params.u = u_m;
            else
                params.u = u_m_noD;
            end
        end
        R = reach(sys,params,options);

        % disturbance
        params.W = W;
        R = reach(sys,params,options);

        % disturbance and sensor noise
        if sys.nrOfOutputs == n
            params.V = V_n;
        elseif sys.nrOfOutputs == y
            params.V = V_y;
        end
        R = reach(sys,params,options);

        % remove fields for clean next iteration
        params = rmfield(params,'U');
        params = rmfield(params,'u');
        params = rmfield(params,'W');
        params = rmfield(params,'V');
        
        if k > 1
            options_cont = rmfield(options_cont,'linAlg');
        end
        
    end
        
end


% all checks ok
res = true;

end


%------------- END OF CODE -------------
