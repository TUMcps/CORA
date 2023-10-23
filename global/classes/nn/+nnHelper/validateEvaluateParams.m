function evParams = validateEvaluateParams(evParams)
% validateEvaluateParams - checks input and sets default values for the
%    evParams struct for the neuralNetwork/evaluate function.
%
% Syntax:
%    evParams = nnHelper.validateEvaluateParams(evParams)
%
% Inputs:
%    evParams - struct (see neuralNetwork/evaluate)
%
% Outputs:
%    evParams - updated evParams
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate

% Authors:       Tobias Ladner
% Written:       29-November-2022
% Last update:   16-February-2023 (poly_method)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
persistent defaultFields
if isempty(defaultFields)
    defaultFields = {
       'bound_approx', true;
       'poly_method', 'regression';
       'num_generators', [];
       'max_gens_post', @(evParams) evParams.num_generators * 100;
       'add_approx_error_to_GI', false;
       'plot_multi_layer_approx_info', false;
       'reuse_bounds', false;
       'max_bounds', 5;
       'do_pre_order_reduction', true;
       'remove_GI', @(evParams) ~evParams.add_approx_error_to_GI;
       'force_approx_lin_at', Inf;
       'propagate_bounds', @(evParams) evParams.reuse_bounds;
       'sort_exponents', false;
       'maxpool_type', 'project';
       'order_reduction_sensitivity', false;
       'graph', graph();
    };
end

% set default value of fields if required
for i=1:size(defaultFields, 1)
    field = defaultFields{i, 1};
    if ~isfield(evParams, field)
        fieldValue = defaultFields{i, 2};
        if isa(fieldValue, "function_handle")
            fieldValue = fieldValue(evParams);
        end
        evParams.(field) = fieldValue;
    end
end

% check fields
if CHECKS_ENABLED
    structName = inputname(1);
    
    % bound_approx
    aux_checkFieldClass(evParams, 'bound_approx', ...
        {'logical', 'string'}, structName);
    if isa(evParams.bound_approx, 'string')
        aux_checkFieldStr(evParams, 'bound_approx', {'sample'}, structName)
        warning("Choosing Bound estimation '%s' does not lead to safe verification!", ...
            evParams.bound_approx);
    end
    
    % poly_method
    validPolyMethod = {'regression', 'ridgeregression', ...
    'throw-catch', 'taylor', 'singh'};
    aux_checkFieldStr(evParams, 'poly_method', validPolyMethod, structName);
    
    % num_generators
    aux_checkFieldClass(evParams, 'num_generators', {'double'}, structName);
    
    % add_approx_error_to_GI
    aux_checkFieldClass(evParams, 'add_approx_error_to_GI', ...
        {'logical'}, structName);
    
    % plot_multi_layer_approx_info
    aux_checkFieldClass(evParams, 'plot_multi_layer_approx_info', ...
        {'logical'}, structName);
    
    % reuse_bounds
    aux_checkFieldClass(evParams, 'reuse_bounds', {'logical'}, structName);
    
    % max_bounds
    aux_checkFieldClass(evParams, 'max_bounds', {'double'}, structName);
    
    % reuse_bounds
    aux_checkFieldClass(evParams, 'do_pre_order_reduction', ...
        {'logical'}, structName);
    
    % max_gens_post
    aux_checkFieldClass(evParams, 'max_gens_post', {'double'}, structName);
    
    % remove_GI
    aux_checkFieldClass(evParams, 'remove_GI', {'logical'}, structName);
    
    % force_approx_lin_at
    aux_checkFieldClass(evParams, 'force_approx_lin_at', {'double'}, structName);
    
    % propagate_bounds
    aux_checkFieldClass(evParams, 'propagate_bounds', {'logical'}, structName);
    
    % sort_exponents
    aux_checkFieldClass(evParams, 'sort_exponents', {'logical'}, structName);
    
    % maxpool_type
    aux_checkFieldStr(evParams, 'maxpool_type', {'project', 'regression'}, structName);
    
    % order_reduction_sensitivity
    aux_checkFieldClass(evParams, 'order_reduction_sensitivity', {'logical'}, structName);
    
    % graph
    aux_checkFieldClass(evParams, 'graph', {'graph'}, structName);
end

end


% Auxiliary functions -----------------------------------------------------


% see also inputArgsCheck, TODO integrate

function aux_checkFieldStr(evParams, field, admissibleValues, structName)
    if CHECKS_ENABLED
        fieldValue = evParams.(field);
        if ~(isa(fieldValue, 'string') || isa(fieldValue, 'char')) || ...
            ~ismember(fieldValue, admissibleValues)
            throw(CORAerror('CORA:wrongFieldValue', ...
                aux_getName(structName, field), admissibleValues))
        end
    end
end

function aux_checkFieldClass(evParams, field, admissibleClasses, structName)
    if CHECKS_ENABLED
        if ~ismember(class(evParams.(field)), admissibleClasses)
            throw(CORAerror('CORA:wrongFieldValue', ...
                aux_getName(structName, field), admissibleClasses))
        end
    end
end

function msg = aux_getName(structName, field)
    msg = sprintf("%s.%s", structName, field);
end

% ------------------------------ END OF CODE ------------------------------
