function mergedBC = ProductMerge(instances)
% ProductMerge - Combines a tree of component instances into a single
%    instance, by building the hybrid automaton product of all
%    Base-Components.
%    note: Function does not check whether input is a tree structure.
%          All contained Base Components are assumed to be tree leaves.
%
% Syntax:
%    mergedBC = ProductMerge(instances)
%
% Inputs:
%    instances - cell array of instantiated components
%        (output of InstantiateComponents)
%
% Outputs:
%    mergedBC - Base-Component instance of the hybrid automaton product
%
% Example:
%    mergedInstance = ProductMerge(instanceArray)
%
% Other m-files required: none
% Subfunctions: automatonProduct
% MAT-files required: none
%
% See also: none

% Authors:       Johann Schoepfer
% Written:       09-April-2018 
% Last update:   09-April-2018 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% output on command window
disp("merge components into flat automaton using automaton product...");

% read out all base components
isBC = cellfun(@(x) ~x.isNetwork,instances,'UniformOutput',true);
BCinstances = instances(isBC);

% Since the automaton product is associative, internal tree structures can
% be ignored: compute (((BC1 X BC2) X BC3) X BC4) ...
mergedBC = BCinstances{1};
for i = 2:numel(BCinstances)
    mergedBC = automatonProduct(mergedBC,BCinstances{i});
end

% add the variable list of the root instance
mergedBC.listOfVar = instances{1}.listOfVar;

% output on command window
fprintf("done. %i BC instances merged into monolithic automaton instance.\n",...
        numel(BCinstances));

% ------------------------------ END OF CODE ------------------------------
