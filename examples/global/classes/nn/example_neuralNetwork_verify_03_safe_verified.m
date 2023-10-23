function res = example_neuralNetwork_verify_03_safe_verified()
% example_neuralNetwork_verify_03_safe_verified - example for 
%    neural network verification using an safe set as specification
%
% Syntax:
%    res = example_neuralNetwork_verify_03_safe_verified()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Authors:       Tobias Ladner
% Written:       01-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng(1)

% get a neural network
nn = neuralNetwork.generateRandom( ...
    'NrInputs', 2, ...
    'NrOutputs', 2, ...
    'ActivationFun', 'tanh', ...
    'NrLayers', 3, ...
    'NrHiddenNeurons', 6 ...
);

% input set
X0 = interval.generateRandom('Dimension', 2);

% specification
spec = specification(interval([-1;-0.5], [1;1]), 'safeSet');

% run verification
doPlot = length(dbstack) == 1; % only plot if called directly.
[isVerified, x] = nn.verify(X0, spec, 'Plot',doPlot,'Verbose',doPlot);

if doPlot
    if isVerified
        disp("Verification was successful!")
    else
        disp("Verification was not succesful. Counter-example:")
        disp("Counterexample:")
        disp(x)
        disp("In output space:")
        y = nn.evaluate(x);
        disp(y)
    end
end

% verification should succeed
res = isVerified;
    
% ------------------------------ END OF CODE ------------------------------
