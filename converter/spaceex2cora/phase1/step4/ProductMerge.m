function [automaton] = ProductMerge(instances)
% ProductMerge - Combines a tree of component instances into a single
%    instance, by building the hybrid automaton product of all
%    Base-Components.
%    (NOTE: Function does not check, wether input is a tree structure.
%       All contained Base-Components are assumed to be tree leaves.)
%
% Syntax:  
%    automaton = ProductMerge(instances)
%
% Inputs:
%    instances - cell array of instantiated components
%        (output of InstantiateComponents)
%
% Outputs:
%    automaton - Base-Component instance of the hybrid automaton product
%
% Example: mergedInstance = ProductMerge(instanceArray)
%
% Other m-files required: none
% Subfunctions: automatonProduct
% MAT-files required: none
%
% See also: none

% Author: Johann Schöpfer
% Written: 09-April-2018 
% Last update: 09-April-2018 
% Last revision: ---

%------------- BEGIN CODE --------------

disp("merge using ProductMerge...");

% first find all Base Components
isBC = false(size(instances));

for i = 1:numel(instances)
    isBC(i) = ~(instances{i}.isNetwork);
end

BCinstances = instances(isBC);

% Since the automaton product is associative,
% internal tree structures can be ingored.
% compute (((BC1 X BC2) X BC3) X BC4) ...
mergedBC = BCinstances{1};

for i = 2:numel(BCinstances)
    mergedBC = automatonProduct(mergedBC,BCinstances{i});
end

% add the variable list of the root instance
mergedBC.listOfVar = instances{1}.listOfVar;

% return result
automaton = mergedBC;

fprintf("done. %i BC instances merged into monolithic automaton instance.\n",...
        numel(BCinstances));

end

%------------- END OF CODE -------------