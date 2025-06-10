function res = test_nn_neuralNetwork_explain()
% test_nn_neuralNetwork_explain - unit test function for 
%     neuralNetwork/explain
%
% Syntax:
%    res = test_nn_neuralNetwork_explain()
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
% See also: -

% Authors:       Tobias Ladner
% Written:       07-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct network
W1 = [ -9 -8 -7 ; 10 -6 0 ; -6 2 5 ; 4 4 -8 ; -5 -8 2 ; 0 6 2 ; -7 10 -2 ; 0 8 6 ; 1 -3 -2 ; 3 9 2 ];
W2 = [ 3 6 -5 3 -6 2 6 2 -4 8 ; 4 1 7 -3 -4 4 2 0 2 -1 ; -3 9 1 5 10 9 1 4 -6 -7 ];
nn = neuralNetwork({nnLinearLayer(W1),nnSigmoidLayer(),nnLinearLayer(W2)});

% construct input
x = [1;2;3];
label = 1;

% compute explanation
verbose = false;
epsilon = 0.2;

% method: standard
method = 'standard';
[idxFreedFeatsStandard] = nn.explain(x,label,epsilon,'InputSize',[3,1,1],'Method',method,'Verbose',verbose);
assert(isequal(idxFreedFeatsStandard,[3,2]));

% method: abstract+refine
method = 'abstract+refine';
[idxFreedFeatsStandard] = nn.explain(x,label,epsilon,'InputSize',[3,1,1],'Method',method,'Verbose',verbose);
assert(isequal(idxFreedFeatsStandard,[3,2]));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
