function BClist = createBClist(instances)
% createBClist - Create list of components running in parallel within a
%    parallel hybrid automaton;
%    note: function does not check whether input is a tree structure.
%          All contained Base-Components are assumed to be tree leaves.
%
% Syntax:
%    BClist = createBClist(instances)
%
% Inputs:
%    instances - cell array of instantiated components
%        (output of InstantiateComponents)
%
% Outputs:
%    BClist - array of BC objects
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       29-November-2018 
% Last update:   --- 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% output on command window
disp("create components of parallel hybrid automaton...");

% read out all base components
isBC = cellfun(@(x) ~x.isNetwork,instances,'UniformOutput',true);
BCinstances = instances(isBC);

% init struct
BClist = struct;

% add information of every base component to struct
for i=1:numel(BCinstances)
    BClist(i).id = BCinstances{i}.id;
    BClist(i).name = BCinstances{i}.name;
    BClist(i).States = BCinstances{i}.States;
    BClist(i).listOfVar = BCinstances{i}.listOfVar;
    BClist(i).listOfLabels = BCinstances{i}.listOfLabels;
    BClist(i).listOfConstants = BCinstances{i}.listOfConstants;
end

% output on command window
fprintf("done. created list containing %i BC instances.\n",...
        numel(BCinstances));

% ------------------------------ END OF CODE ------------------------------
