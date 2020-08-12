function res = test_optionsValidFields(~)
% test_optionsValidFields - test for the correctness of
% - getValidFields.m
% - getDefaultOption.m
% 
% this unit test checks ...
% 1) if valid fields from options check do not intersect with each other
% 2) if options declared as default all have default value
%
% Syntax:  
%    res = test_optionsValidFields
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: -

% Author:       Mark Wetzlinger
% Written:      12-Aug-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% different kinds of options
types = {'mand';'opt';'def'};
% ...assumption: mand always non-empty

% comparisons:
comps = [1 2; 1 3; 2 3];

% different systems
allsys = {'contDyn_sim';'linearSys';'linearSysDT';'linParamSys';'nonlinDASys';'nonlinearSysDT';...
    'nonlinearSys';'nonlinParamSys';'hybridAutomaton_reach';'hybridAutomaton_sim';...
    'parallelHybridAutomaton_reach';'parallelHybridAutomaton_sim'};

for s=1:length(allsys)
    
    % choose system for which options are to be checked
    sys = allsys{s};
    
    % get groups of valid fields
    for i=1:length(types)
        try
            validFields.(types{i}) = getValidFields(sys,types{i});
        catch
            res = false;
            error("defined system class does not exist");
        end
    end

    % check groups against each other
    % error if same option in different groups (mand, opt, def)
    for c=1:length(comps)
        group1 = validFields.(types{comps(c,1)});
        group2 = validFields.(types{comps(c,2)});
        if ~isempty(group1) && ~isempty(group2)
            if any(cellfun(@any, ...
                cellfun(@(fields) strcmp(fields,group2), group1, 'UniformOutput', false)))
                res = false;
                error(sys + ": validFields of '" + types(comps(c,1)) + "' and '" + ...
                    types(comps(c,2)) + "' intersect");
            end
        end
    end
    
    % check if all options declared as default really have a default value
    if ~isempty(validFields.def)
        try
            cellfun(@getDefaultOption, validFields.def, 'UniformOutput', false);
        catch
            error(sys + ": some default option does not have a default value");
        end
    end
    
    disp("options for system class '" + sys + "' ok");

end

res = true;


end

%------------- END OF CODE --------------
