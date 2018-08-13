function [Rguard,activeGuards,minIndex,maxIndex] = guardIntersect(obj,guards,setIndices,R,options)
% detGuardIntersect - intersects the reachable sets with potential guard sets
% and returns enclosing zonotopes for each guard set for the deterministic
% case
%
% Syntax:  
%    ---
%
% Inputs:
%    obj - location object
%    R - cell array of reachable sets
%
% Outputs:
%    Z - cell array of enclosing zonotopes for each transition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      08-May-2007 
% Last update:  21-September-2007
%               26-March-2008
%               11-October-2008
%               24-July-2009
%               20-October-2010
%               30-July-2016
%               23-November-2017
%               20-April-2018 (intersect guard sets with invariant)
% Last revision:---

%------------- BEGIN CODE --------------

%initialize values
Rguard=[]; activeGuards=[]; minIndex=[]; maxIndex=[];
relevantIndex=unique(setIndices);

%check if guard sets exist
if ~isempty(obj.transition) && ~isempty(relevantIndex)
    
    %initialize intersections
    N = length(obj.transition);
    
    % determine guards that got hit
    guardInd = unique(guards);
    
    % extract guard sets
    Pguard = cell(N,1);

    for i=1:length(guardInd)
        Pguard{guardInd(i)} = get(obj.transition{guardInd(i)},'guard');
    end
    
    
    
    % selected method for the calculation of the intersection
    if strcmp(options.guardIntersect,'polytope')

        %build polytopes of reachable sets    
        for iSet=1:length(relevantIndex)
            if ~iscell(R{iSet}) % no split sets
                [Pset{relevantIndex(iSet)}]=enclosingPolytope(R{relevantIndex(iSet)},options);
            else % sets are split
                for iSubSet=1:length(R{iSet})
                    [Pset{relevantIndex(iSet)}{iSubSet}]=enclosingPolytope(R{relevantIndex(iSet)}{iSubSet},options);
                end
            end
        end
        
        % intersection of candidates
        funcAnd = @(x,y) x & y;
        funcIsempty = @(x) isempty(x);
        [minIndex,maxIndex,Pint,activeGuards] = calcIntersection(Pset,Pguard,guards,setIndices,funcAnd,funcIsempty);
        
        % compute enclosing zonotopes
        Rguard = cell(length(Pint),1);
        
        for i=1:length(activeGuards)
            Rguard{i}{1} = enclosePolytopes(obj,Pint{i}, options);
        end 


        
        
    elseif strcmp(options.guardIntersect,'zonoGirard')
        
        % group the reachable sets which intersect guards
        [minIndex,maxIndex,P,activeGuards] = groupSets(R,guards,setIndices);
        
        % determine multiple suitable orthogonal basis for each group
        D = cell(length(activeGuards),1);

        for i = 1:length(activeGuards)
            D{i} = calcBasis(obj,P{i},Pguard{activeGuards(i)},options);
        end
        
        % calculate the intersection
        Pint = cell(size(P));
        
        for i = 1:length(P)
           
           Pint{i} = cell(length(P{i}),1);
            
           for j = 1:length(P{i})
              if ~iscell(P{i}{j})       % no split sets
                  Pint{i}{j} = intersectHyperplaneGirard(P{i}{j},Pguard{activeGuards(i)},D{i});
              else                      % sets are split
                  for k = 1:length(P{i}{j})              
                      Pint{i}{j}{k} = intersectHyperplaneGirard(P{i}{j}{k},Pguard{activeGuards(i)},D{i});
                  end
                  
                  Pint{i}{j} = Pint{i}{j}(~cellfun('isempty',Pint{i}{j})); 
              end
           end 
        end
        
        % remove all empty intersections
        [Pint,D,minIndex,maxIndex,activeGuards] = removeEmptySets(Pint,D,minIndex,maxIndex,activeGuards);
        
        % compute enclosing zonotopes
        Rguard = cell(length(Pint),1);
        
        for i=1:length(activeGuards)
            Rguard{i}{1} = encloseAlignedIntervals(Pint{i},D{i},Pguard{activeGuards(i)});
        end 
    else
        error('Wrong value for setting options.guardIntersect!');  
    end
end
end



% Auxiliary Functions -----------------------------------------------------

function [minInd,maxInd,Pint,guards] = calcIntersection(Pset,Pguard,guards,setIndices,funcAnd,funcIsempty)
% calculate the intersections between guard sets and reachable sets

    setIndicesGuards = cell(max(guards),1);
    Pint = cell(max(guards),1);
    
    % loop over all sets which hitted a guard set
    for i=1:length(guards)
        
        %intersect polytopes
        if ~iscell(Pset{setIndices(i)}) % no split sets
            intersection = funcAnd(Pset{setIndices(i)},Pguard{guards(i)});
            %intersection empty?
            if ~funcIsempty(intersection)
                setIndicesGuards{guards(i)} = [setIndicesGuards{guards(i)};setIndices(i)];
                Pint{guards(i)}{end+1}=intersection;           
            end
        else % sets are split
            for iSubSet=1:length(Pset{setIndices(i)})
                intersection{iSubSet} = funcAnd(Pset{setIndices(i)}{iSubSet},Pguard{guards(i)});
                %intersection empty?
                if ~funcIsempty(intersection{iSubSet})
                    notEmpty(iSubSet) = 1;
                else
                    notEmpty(iSubSet) = 0;
                end
            end
            % one intersection was not empty
            if any(notEmpty)
                setIndicesGuards{guards(i)} = [setIndicesGuards{guards(i)},setIndices(i)];
                Pint{guards(i)}{end+1}=[];
                for iSubSet=1:length(Pset{setIndices(i)})
                    Pint{guards(i)}{end}{iSubSet}=intersection{iSubSet}; 
                end            
            end
        end
    end
    
    % remove empty entries (guard sets that have not been hit)
    guards = unique(guards);
    setIndicesGuards = setIndicesGuards(guards);
    Pint = Pint(guards);
    
    % remove all empty guard instersections (approximate collision test 
    % reported intersection, but exact collision test showed that no
    % intersection occurs)
    ind = find(~cellfun(@isempty,setIndicesGuards));
    setIndicesGuards = setIndicesGuards(ind);
    guards = guards(ind);
    Pint = Pint(ind);
    
    % Remove gaps in the set-index vector of each intersection
    [minInd,maxInd,Pint,guards] = removeGaps(setIndicesGuards,guards,Pint);
    
end

function [minInd,maxInd,P,guards] = groupSets(Pset,guards,setIndices)
% group the reachable sets which intersect guard sets. The sets in one
% group all intersect the same guard set and are located next to each other

    % initialization
    guardInd = unique(guards);
    setIndicesGuards = cell(length(guardInd),1);
    P = cell(length(guardInd),1);
    
    % Step 1: group accoring to hitted guard sets 
    for i = 1:length(guardInd)
        ind = find(guards == guardInd(i));
        setIndicesGuards{i} = setIndices(ind);
        P{i} = Pset(setIndices(ind));
    end
    
    % Step 2: group accoring to the location (neigbouring sets together)
    [minInd,maxInd,P,guards] = removeGaps(setIndicesGuards,guardInd,P);

end

function [minInd,maxInd,Pint,guards] = removeGaps(setIndicesGuards,guards,Pint)
% Remove gaps in the set-index vector of each intersection

    % split all guard intersections with gaps between the set-indices into
    % multiple different intersections (needed if one guard set is hit
    % multiple times at different points in time)
    counter = 1;
    
    while counter <= length(guards)
       
        setIndices = setIndicesGuards{counter};
        
        for i = 1:length(setIndices)-1
            
            % check if a gap occurs 
            if setIndices(i+1) ~= setIndices(i)+1
               % add first part of the intersection (=gap free) to the 
               % beginning of the list
               setIndicesGuards = [{setIndices(1:i)};setIndicesGuards];
               Pint = [{Pint{counter}(1:i)};Pint];
               guards = [guards(counter);guards];
               
               % add second part of the intesection (possibly contains
               % further gaps) to the part of the list that is not finished
               % yet
               setIndicesGuards{counter+1} = setIndices(i+1:end);
               Pint{counter+1} = Pint{counter+1}(i+1:end);
               
               break;
            end
        end   
        
        counter = counter + 1;
    end
    
    % determine minimum and maximum set-index for each intersection
    minInd = cellfun(@(x) x(1),setIndicesGuards);
    maxInd = cellfun(@(x) x(end),setIndicesGuards);
    
end

function [Pint,D,minIndex,maxIndex,activeGuards] = removeEmptySets(Pint,D,minIndex,maxIndex,activeGuards)
% remove all sets for which the intersection with the guard set turned out
% to be empty

    counter = 1;
    
    % loop over all groups of intersecting sets
    while counter <= length(Pint)
    
        % loop over all single sets that belong to the group
        for i = 1:length(Pint{counter})
            if isempty(Pint{counter}{i})
                
                % split groups if there are further sets in the group
                if i < length(Pint{counter})
                    D{end+1} = D{counter};
                    activeGuards = [activeGuards;activeGuards(counter)];
                    Pint{end+1} = Pint{counter}(i+1:end);
                    minIndex = [minIndex;minIndex(counter)+i];
                    maxIndex = [maxIndex;maxIndex(counter)];
                end
                
                % create a new group from the previous sets if the removed
                % set is not the first one
                if i ~= 1
                   Pint{counter} = Pint{counter}(1:i-1);
                   maxIndex(counter) = minIndex(counter) + i -2;                
                else
                   Pint{counter} = [];
                end
                
                break;
            end
        end
        
        counter = counter + 1;      
    end
    
    % remove all empty sets
    ind = ~cellfun('isempty',Pint);
    
    Pint = Pint(ind);
    D = D(ind);
    minIndex = minIndex(ind);
    maxIndex = maxIndex(ind);
    activeGuards = activeGuards(ind);
    
end


%------------- END OF CODE --------------