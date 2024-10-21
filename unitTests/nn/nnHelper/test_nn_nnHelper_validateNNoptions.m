function res = test_nn_nnHelper_validateNNoptions()
% test_nn_nnHelper_validateNNoptions - tests the params validation for
%    the set-based neural network validation.
%
% Syntax:
%    res = test_nn_nnHelper_validateNNoptions()
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
% See also: nnHelper/validateNNoptions

% Authors:       Tobias Ladner
% Written:       29-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enter empty
options = struct;
options = nnHelper.validateNNoptions(options);

% enter full
options = nnHelper.validateNNoptions(options);

% test fields
options = struct; % test default
options = nnHelper.validateNNoptions(options); 
assert(strcmp(options.nn.poly_method, 'regression'));

options.nn.bound_approx = false;
options = nnHelper.validateNNoptions(options);
assert(options.nn.bound_approx == false);

options.nn.num_generators = 1000;
options = nnHelper.validateNNoptions(options);
assert(options.nn.num_generators == 1000);

options.nn.max_gens_post = 100;
options = nnHelper.validateNNoptions(options);
assert(options.nn.max_gens_post == 100);

options.nn.add_approx_error_to_GI = true;
options = nnHelper.validateNNoptions(options);
assert(options.nn.add_approx_error_to_GI == true);

options.nn.plot_multi_layer_approx_info = true;
options = nnHelper.validateNNoptions(options);
assert(options.nn.plot_multi_layer_approx_info == true);

options.nn.poly_method = 'ridgeregression';
options = nnHelper.validateNNoptions(options);
assert(strcmp(options.nn.poly_method, 'ridgeregression'));

options.nn.reuse_bounds = true;
options = nnHelper.validateNNoptions(options);
assert(options.nn.reuse_bounds == true);

options.nn.max_bounds = 2;
options = nnHelper.validateNNoptions(options);
assert(options.nn.max_bounds == 2);

options.nn.do_pre_order_reduction = false;
options = nnHelper.validateNNoptions(options);
assert(options.nn.do_pre_order_reduction == false);

options.nn.remove_GI = false;
options = nnHelper.validateNNoptions(options);
assert(options.nn.remove_GI == false);

options.nn.force_approx_lin_at = 2;
options = nnHelper.validateNNoptions(options);
assert(options.nn.force_approx_lin_at == 2);

% test error message
options = struct;
options.nn.poly_method = 'unknown';
assertThrowsAs(@nnHelper.validateNNoptions,'CORA:wrongFieldValue',options);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
