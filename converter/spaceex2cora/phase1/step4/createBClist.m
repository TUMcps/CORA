function [BClist] = createBClist(instances)
% createBClist - ?
%    (NOTE: Function does not check, whether input is a tree structure.
%       All contained Base-Components are assumed to be tree leaves.)
%
% Syntax:  
%    HAlist = createBClist(instances)
%
% Inputs:
%    instances - cell array of instantiated components
%        (output of InstantiateComponents)
%
% Outputs:
%    HAlist - array of BC objects
%
% Other m-files required: none
% Subfunctions: ADD IF NECESSARY
% MAT-files required: none
%
% See also: none

% Author: Mark Wetzlinger
% Written: 29-November-2018 
% Last update: --- 
% Last revision: ---

%------------- BEGIN CODE --------------

disp("create components of parallel hybrid automaton using createBClist...");

% first find all BaseComponents
isBC = false(size(instances));

for i = 1:numel(instances)
    isBC(i) = ~(instances{i}.isNetwork);
end

BCinstances = instances(isBC);

BClist = struct;

% create instance for every BaseComponent
for i=1:numel(BCinstances)
    BClist(i).id = BCinstances{i}.id;
    BClist(i).name = BCinstances{i}.name;
    BClist(i).States = BCinstances{i}.States;
    BClist(i).listOfVar = BCinstances{i}.listOfVar;
end

fprintf("done. created list containing %i BC instances.\n",...
        numel(BCinstances));

end

%------------- END OF CODE -------------
