function res = example_neuralNetwork_verify_01_unsafe_verified()
% example_neuralNetwork_verify_01_unsafe_verified - example for 
%    neural network verification using an unsafe set as specification
%
% Syntax:
%    res = example_neuralNetwork_verify_01_unsafe_verified()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-December-2020
% Last update:   30-November-2022 (removed neuralNetworkOld)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng(2)

% get a neural network
nn = neuralNetwork.generateRandom( ...
    'NrInputs', 2, ...
    'NrOutputs', 2, ...
    'ActivationFun', 'tanh', ...
    'NrLayers', 3 ...
);

% input set
X0 = interval.generateRandom('Dimension', 2);

% specification
spec = specification(interval([0;-0.25], [1;0.65]), 'unsafeSet');

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
