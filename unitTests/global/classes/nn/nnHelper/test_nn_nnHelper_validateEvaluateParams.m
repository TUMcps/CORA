function res = test_nn_nnHelper_validateEvaluateParams()
% test_nn_nnHelper_validateEvaluateParams - tests the params validation for
%    the set-based neural network validation.
%
% Syntax:
%    res = test_nn_nnHelper_validateEvaluateParams()
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
% See also: nnHelper/validateEvaluateParams

% Authors:       Tobias Ladner
% Written:       29-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% enter empty
evParams = struct;
evParams = nnHelper.validateEvaluateParams(evParams);

% enter full
evParams = nnHelper.validateEvaluateParams(evParams);

% test fields
evParams = struct; % test default
evParams = nnHelper.validateEvaluateParams(evParams); 
res = res && strcmp(evParams.poly_method, 'regression');

evParams.bound_approx = false;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.bound_approx == false;

evParams.num_generators = 1000;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.num_generators == 1000;

evParams.max_gens_post = 100;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.max_gens_post == 100;

evParams.add_approx_error_to_GI = true;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.add_approx_error_to_GI == true;

evParams.plot_multi_layer_approx_info = true;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.plot_multi_layer_approx_info == true;

evParams.poly_method = 'ridgeregression';
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && strcmp(evParams.poly_method, 'ridgeregression');

evParams.reuse_bounds = true;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.reuse_bounds == true;

evParams.max_bounds = 2;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.max_bounds == 2;

evParams.do_pre_order_reduction = false;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.do_pre_order_reduction == false;

evParams.remove_GI = false;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.remove_GI == false;

evParams.force_approx_lin_at = 2;
evParams = nnHelper.validateEvaluateParams(evParams);
res = res && evParams.force_approx_lin_at == 2;


% test error message
try
    evParams = struct;
    evParams.poly_method = 'unknown';
    evParams = nnHelper.validateEvaluateParams(evParams);
    % should have trown an error
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
