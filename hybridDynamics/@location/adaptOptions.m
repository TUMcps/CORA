function [params,options] = adaptOptions(obj,options)
% adaptOptions - extract the settings for continous reachability from the
%                settings for hybrid reachability
%
% Syntax:  
%    [params,options] = adaptOptions(obj,options)
%
% Inputs:
%    obj - location object
%    options - settings for hybrid reachability analysis
%
% Outputs:
%    params - parameter for continuous system reachability analysis
%    options - settings for continuous system reachability analysis
%
% See also: location/reach

% Author:       Niklas Kochdumper, Mark wetzlinger
% Written:      09-June-2020 
% Last update:  08-February-2021 (MW, adapt to new params/options check)
% Last revision: ---

%------------- BEGIN CODE --------------

% read params and options from contDynamics
[paramsList,optionsList] = feval(['config_' class(obj.contDynamics) '_reach'],...
    obj.contDynamics,options,options);

% copy parameter from options
params.U = options.U;
params.tFinal = options.tFinal;
if isfield(options,'paramInt')
    params.paramInt = options.paramInt; 
    options = rmField(options,'paramInt');
end

% remove parameter from options
for i=1:length(paramsList.name)
    if isfield(options,paramsList.name{i})
        options = rmfield(options,paramsList.name{i}); 
    end
end 

% remove all options not for contDynamics class
optFields = fields(options);
for i=1:length(optFields)
    if ~ismember(optionsList.name,optFields{i})
        options = rmfield(options,optFields{i}); 
    end
end

% % remove parameters from options
% options = rmfield(options,'R0');
% options = rmfield(options,'tFinal');
% options = rmfield(options,'tStart');
% options = rmfield(options,'startLoc');
% options = rmfield(options,'U');
% 
% % remove redundant fields
% valFields = getValidFields(class(obj.contDynamics));
% list = fields(options);
% 
% for i = 1:length(list)
%     if ~ismember(list{i},valFields)
%         options = rmfield(options,list{i}); 
%     end
% end
    
% compute time point solution for linear parametric systems
if isa(obj.contDynamics,'linParamSys')
    options.compTimePoint = true;
end

end

%------------- END OF CODE --------------