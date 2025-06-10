function res = example_neuralNetwork_verify_safe()
% example_neuralNetwork_verify_safe - example for the verification of a 
%    neural networks using the function neuralNetwork/verify.
%
% Syntax:
%    res = example_neuralNetwork_verify_safe()
%
% Inputs:
%    -
%
% Outputs:
%    res - string, verification result 
%       ['VERIFIED','COUNTEREXAMPLE','UNKNOWN']
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       18-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng('default')

verbose = true;
% Specify model and specification path.
modelPath = 'ACASXU_run2a_1_2_batch_2000.onnx';
specPath = 'prop_1.vnnlib';
timeout = 2;
% Load model and specification.
[nn,x,r,A,b,safeSet,options] = aux_readModelAndSpecs(modelPath,specPath);
% Do verification.
timerVal = tic;
[res,x_,y_] = nn.verify(x,r,A,b,safeSet,options,timeout,verbose);
% Print result.
if verbose
    % Print result.
    fprintf('%s -- %s: %s\n',modelPath,specPath,res);
    time = toc(timerVal);
    fprintf('--- Verification time: %.4f / %.4f [s]\n',time,timeout);
end
% Write results.
fprintf(['Result:' newline]);
aux_writeResults(res,x_,y_);

end


% Auxiliary functions -----------------------------------------------------

function [nn,x,r,A,b,safeSet,options] = aux_readModelAndSpecs(modelPath,specPath)
    % Load the model.
    nn = neuralNetwork.readONNXNetwork(modelPath,false,'BSSC');
    % Load specification.
    [X0,specs] = vnnlib2cora(specPath);
    % Compute center and radius of the input set.
    x = 1/2*(X0{1}.sup + X0{1}.inf);
    r = 1/2*(X0{1}.sup - X0{1}.inf);

    % Extract specification.
    if isa(specs.set,'halfspace')
        A = specs.set.c';
        b = -specs.set.d;
    else
        A = specs.set.A;
        b = -specs.set.b;
    end
    safeSet = strcmp(specs.type,'safeSet');

    % Create evaluation options.
    options.nn = struct(...
        'use_approx_error',true,...
        'poly_method','bounds',...'bounds','singh'
        'train',struct(...
            'backprop',false,...
            'mini_batch_size',512 ...
        ) ...
    );
    % Set default training parameters
    options = nnHelper.validateNNoptions(options,true);
    options.nn.interval_center = false;
end

function aux_writeResults(res,x_,y_)
    % Write results.
    if strcmp(res,'VERIFIED')
        % Write content.
        fprintf(['unsat' newline]);
    elseif strcmp(res,'COUNTEREXAMPLE')
        % Write content.
        fprintf(['sat' newline '(']);
        % Write input values.
        for j=1:size(x_,1)
            fprintf(['(X_%d %f)' newline],j-1,x_(j));
        end
        % Write output values.
        for j=1:size(y_,1)
            fprintf(['(Y_%d %f)' newline],j-1,y_(j));
        end
        fprintf(')');
    else
        fprintf(['unknown' newline]);
    end

end

% ------------------------------ END OF CODE ------------------------------
