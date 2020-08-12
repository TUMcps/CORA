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
%    options - settings fro continuous system reachability analysis
%
% See also: location/reach

% Author:       Niklas Kochdumper
% Written:      09-June-2020 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % copy parameter from options
    params.U = options.U;
    params.tFinal = options.tFinal;
    
    if isfield(options,'paramInt')
       params.paramInt = options.paramInt; 
       options = rmField(options,'pararmInt');
    end
    
    % remove parameter from options
    options = rmfield(options,'tFinal');
    options = rmfield(options,'tStart');
    options = rmfield(options,'startLoc');
    options = rmfield(options,'U');
    
    % remove redundant fields
    valFields = getValidFields(class(obj.contDynamics));
    list = fields(options);
    
    for i = 1:length(list)
       if ~ismember(list{i},valFields)
          options = rmfield(options,list{i}); 
       end
    end
    
    % compute time point solution for linear parametric systems
    if isa(obj.contDynamics,'linParamSys')
        options.compTimePoint = 1;
    end
end

%------------- END OF CODE --------------