function [params,options] = adaptOptions(loc,options)
% adaptOptions - extract the settings for continous reachability from the
%    settings for hybrid reachability
%
% Syntax:
%    [params,options] = adaptOptions(loc,options)
%
% Inputs:
%    loc - location object
%    options - settings for hybrid reachability analysis
%
% Outputs:
%    params - parameters for continuous system reachability analysis
%    options - settings for continuous system reachability analysis
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/reach

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       09-June-2020 
% Last update:   08-February-2021 (MW, adapt to new params/options check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read params and options from contDynamics
[paramsList,optionsList] = feval(['config_' class(loc.contDynamics) '_reach'],...
    loc.contDynamics,options,options);

% copy parameter from options
params.U = options.U;
params.tFinal = options.tFinal;
if isfield(options,'paramInt')
    params.paramInt = options.paramInt; 
    options = rmfield(options,'paramInt');
end

% remove model parameters from options
for i=1:length(paramsList)
    if isfield(options,paramsList(i).name)
        options = rmfield(options,paramsList(i).name); 
    end
end 

% read options for contDynamics class
optFields = fields(options);
% we need some additional handling if a field is a struct
addFields = {}; removeIdx = false(length(optFields),1);
for i=1:length(optFields)
    % check if field is a struct, requires one depth more
    if isstruct(options.(optFields{i}))
        fieldsStruct = fields(options.(optFields{i}));
        % loop over all fields one level deeper and add them to the list
        for j=1:length(fieldsStruct)
            addFields = [addFields; {[optFields{i} '.' fieldsStruct{j}]}];
            % remove high-level struct 
            removeIdx(i) = true;
        end
    end
end
optFields = optFields(~removeIdx);
optFields = [optFields; addFields];

% remove all options not for contDynamics class
for i=1:length(optFields)
    if ~ismember({optionsList.name}',optFields{i})
        options = rmfield(options,optFields{i}); 
    end
end
    
% compute time point solution for linear parametric systems
if isa(loc.contDynamics,'linParamSys')
    options.compTimePoint = true;
end

% ------------------------------ END OF CODE ------------------------------
