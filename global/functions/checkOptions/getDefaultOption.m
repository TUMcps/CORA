function defValue = getDefaultOption(option)
% getDefaultOption - return default value for given option
% ...currently only for linearSys
%
% Syntax:
%    defValue = getDefaultOption(option)
%
% Inputs:
%    option   - string <opt> in options.<opt>
%
% Outputs:
%    defValue - default value of given option
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkOptionsReach.m, checkOptionsSimulate.m
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      08-Aug-2019
% Last update:  05-May-2020 (addition of reductionInterval)
% Last revision:---

%------------- BEGIN CODE --------------

% {'name of option #1', default value #1;
%    'name of option #2', default value #2;
%    ..., ...};
optValues = {'tStart', 0;
             'reductionTechnique', 'girard';
             'linAlg', 'standard';
             'verbose', false;
             'reductionInterval', Inf;
             'maxError', Inf; 
             'compTimePoint', true};
         
[row,~] = find(strcmp(optValues,option));
if ~isempty(row)
    defValue = optValues{row,2};
else
    error("There is no default value for params/options." + option + "!");
end

end

%------------- END OF CODE --------------
