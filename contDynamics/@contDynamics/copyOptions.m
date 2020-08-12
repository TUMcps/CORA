function newOptions = copyOptions(obj,options)
% copyOptions - copies options of system
%
% Syntax:  
%    copyOptions(obj,options)
%
% Inputs:
%    obj     - object
%    options - options
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       ---
% Written:      ---
% Last update:  08-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

if isa(obj,'linearSys')
    vf = getValidFields(class(obj),'all',options.linAlg);
else
    vf = getValidFields(class(obj));
end

defaults = getValidFields(class(obj),'def');
defValues = cellfun(@(opt) getDefaultOption(opt), defaults, 'UniformOutput', false);
for i=1:numel(vf)
    if isfield(options,vf{i})
        defIdx = find(strcmp(vf{i},defaults),1);
        if ~isempty(defIdx)
            % is default option -> check for default value
            if ischar(options.(vf{i}))
                if ~strcmp(options.(vf{i}),defValues{defIdx})
                    newOptions.(vf{i}) = options.(vf{i});
                end
            elseif isnumeric(options.(vf{i})) || islogical(options.(vf{i}))
                if ~(options.(vf{i}) == defValues{defIdx})
                    newOptions.(vf{i}) = options.(vf{i});
                end
            else
                warning("Unknown value comparison");
                newOptions.(vf{i}) = options.(vf{i});
            end
        else
            newOptions.(vf{i}) = options.(vf{i});
        end
    end     
end

%------------- END OF CODE --------------