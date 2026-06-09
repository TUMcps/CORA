function res = test_nn_nnActivationLayer_getMergeBuckets()
% test_nn_nnActivationLayer_getMergeBuckets - tests the suitable points for
%    merging neurons per activation
%
% Syntax:
%    res = test_nn_nnActivationLayer_getMergeBuckets()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       23-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

layer = nnSigmoidLayer();
buckets = layer.getMergeBuckets();
res = res && all(buckets == [0, 1]);

layer = nnTanhLayer();
buckets = layer.getMergeBuckets();
res = res && all(buckets == [-1, 1]);

layer = nnReLULayer();
buckets = layer.getMergeBuckets();
res = res && all(buckets == 0);

layer = nnLeakyReLULayer();
buckets = layer.getMergeBuckets();
res = res && all(buckets == []);

% ------------------------------ END OF CODE ------------------------------
