function sys = nonlinearSysDT(obj)
% nonlinearSysDT - transform a neural netowrk to a nonlinarSysDT object
%
% Syntax:
%    sys = nonlinearSysDT(obj)
%
% Inputs:
%    obj - object of class neuralNetwork
%
% Outputs:
%    sys - object of class nonlinearSysDT
%
% References:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Laura Luetzow
% Written:       07-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
assertNarginConstructor([0 1],nargin);

% create function handle for the output of the neural network
f = "x";
dim_u = 0;
for i=1:length(obj.layers)
    if isa(obj.layers{i}, 'nnLinearLayer')
        f = sprintf("net.layers{%d}.W*(%s) + net.layers{%d}.b + u(%d:%d)", i, f, i, dim_u+1, dim_u+size(obj.layers{i}.b,1));
        dim_u = dim_u+size(obj.layers{i}.b,1);
    elseif isa(obj.layers{i}, 'nnElementwiseAffineLayer')
        f = sprintf("(%s) + net.layers{%d}.offset", f, i);
    elseif isa(obj.layers{i}, 'nnTanhLayer')
        f = sprintf("tanh(%s)", f);
    elseif isa(obj.layers{i}, 'nnReLULayer')
        f = sprintf("extractdata(relu(dlarray(%s)))", f);
    elseif isa(obj.layers{i}, 'nnSigmoidLayer')
        f = sprintf("extractdata(sigmoid(dlarray(%s)))", f);
    elseif isa(obj.layers{i}, 'nnSoftmaxLayer')
        % don't add softmax to network function
    else
        throw(CORAerror('CORA:notSupported', 'Layer type is not implemented for conversion'));
    end
end
dim_x = obj.layers{1}.inputSize(1);
dim_y = obj.neurons_out;
func = eval(sprintf("@(x,u,net) %s", f)); % func(x,0) is equal to evaluate(net, x)

% create nonlinearSysDT object
sys = nonlinearSysDT('nnSys',@(x,u) x,1,dim_x,dim_u,@(x,u) func(x,u,obj),dim_y,false);

% set derivatives
deriv_func = @(x,u) aux_computeJacobianNN(sys,x,u); 
sys = sys.setJacobian(deriv_func);
deriv_func_out = @(x,u) aux_computeOutJacobianNN(obj,x,u);
sys = sys.setOutJacobian(deriv_func_out);
end


% Auxiliary functions -----------------------------------------------------

function [A,B] = aux_computeJacobianNN(sys,x,u)
% compute Jacobian for the state, which is assumed to be constant
A = eye(sys.nrOfDims);
B = zeros(sys.nrOfDims,sys.nrOfInputs);
end

function [C,D] = aux_computeOutJacobianNN(net,x,u)
% compute Jacobian of the output
C = calcSensitivity(net,x);
D = [];

if isa(net.layers{end}, 'nnSoftmaxLayer')
    numLayers = length(net.layers)-2;
else
    numLayers = length(net.layers)-1;
end
for i = 1:numLayers
    if isa(net.layers{i}, 'nnLinearLayer')
        if net.layers{i+1}.sensitivity == 1
            D = [D eye(net.neurons_out)];
        else
            D = [D net.layers{i+1}.sensitivity];
        end
    end
end
D = [D eye(net.neurons_out)];
end

% ------------------------------ END OF CODE ------------------------------
