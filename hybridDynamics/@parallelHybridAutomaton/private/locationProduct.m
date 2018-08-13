function res = locationProduct(obj, loc)
% locationProduct - Construct a overall location object from the active 
%                   loctions of the subcomponts with a local automaton
%                   product
%
% Syntax:  
%    res = locationProduct(obj, loc)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    loc - ID's for the active location of each subcomponent
%
% Outputs:
%    res - constructed location object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schöpfer, Niklas Kochdumper
% Written:      08-June-2018  
% Last update:  09-July-2018
% Last revision: ---

%------------- BEGIN CODE --------------

    numComp = length(obj.components);

    % get active location objects
    locList = cell(1,numComp);

    for i = 1:numComp
        locList{i} = getLocation(obj.components{i},loc{i});
    end


    % preallocate arrays
    locNames = cell(1,numComp);
    invList = cell(1,numComp);
    flowList = cell(1,numComp);
    transList = cell(1,numComp);

    % get properties from all active components
    for iLocation = 1:numComp
        locNames{iLocation} = get(locList{iLocation},'name');   
        invList{iLocation} = get(locList{iLocation},'invariant');
        flowList{iLocation} = get(locList{iLocation},'contDynamics');
        transList{iLocation} = get(locList{iLocation},'transition');
    end

    % create merged name in tuple notation: '(s1,s2,s3)'
    joinedName = join(locNames,',');
    mergedLocName = ['(' joinedName{1} ')'];

    % merge invariants
    mergedInvSet = mergeInvariant(obj, invList);

    % merge flows (= continious dynamics)
    mergedFlow = mergeFlows(obj, flowList);

    % merge transitions
    mergedTransSets = mergeTransitionSets(obj, transList, loc, mergedInvSet);

    % construct resulting location object
    res = location(mergedLocName,loc,mergedInvSet,mergedTransSets,mergedFlow);

end

%------------- END OF CODE --------------