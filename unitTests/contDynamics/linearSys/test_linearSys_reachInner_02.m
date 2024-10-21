function res = test_linearSys_reachInner_02
% test_linearSys_reachInner_02 - compare the computed inner approximation
%    with a computed outer approximation in different cases for dynamics
%    and uncertain inputs
%
% Syntax:
%    res = test_linearSys_reachInner_02
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Mark Wetzlinger
% Written:       06-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dynamics
A = [-1 -4; 4 -1]; B = [1; 0.5]; c = [-0.2; 0.5]; C = [1 -1];
sys{1} = linearSys('2D',A,B);
sys{2} = linearSys('2D',A,B,c);
sys{3} = linearSys('2D',A,B,[],C);

% input sets
U{1} = zonotope(0);
U{2} = zonotope(0.2);
U{3} = zonotope(0.5,0.2);

for i=3:length(sys)
    for j=1:length(U)
        [Router,Rinner] = aux_computeOuterInnerApprox(sys{i},U{j});
        for k=1:length(Rinner.timePoint.set)
            % requirement: same time step size for inner and outer approximation
            assertLoop(contains(Router.timePoint.set{k},Rinner.timePoint.set{k}),i,j,k);
        end
    end
end

% test completed
res = true;

end 


% Auxiliary functions -----------------------------------------------------

function [Router,Rinner] = aux_computeOuterInnerApprox(sys,U)
% set model parameters and reachability options here, then run computation
% of outer approximation and inner approximation

    % model parameters
    params.tFinal = 1;
    params.R0 = zonotope([1;1],0.1*eye(2));
    params.U = U;
    
    % reachability options
    optionsInner.timeStep = 0.01;
    optionsInner.zonotopeOrder = 20;
    optionsOuter.timeStep = 0.01;
    optionsOuter.zonotopeOrder = 20;
    optionsOuter.taylorTerms = 5;
    optionsOuter.linAlg = 'standard';
    
    % compute outer and inner approximation
    Router = reach(sys,params,optionsOuter);
    Rinner = reachInner(sys,params,optionsInner);

end

% ------------------------------ END OF CODE ------------------------------
