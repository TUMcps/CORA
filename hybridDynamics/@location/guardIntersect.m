function [Rguard,actGuards,minInd,maxInd] = ...
    guardIntersect(loc,guards,setInd,setType,Rcont,options)
% guardIntersect - computes an enclosure of the intersection between the
%    reachable set and the guard sets
%
% Syntax:
%    [Rguard,actGuards,minInd,maxInd] = 
%       guardIntersect(loc,guards,setInd,setType,Rcont,options)
%
% Inputs:
%    loc - location object
%    guards - cell array containing the guard sets that have been hit
%    setInd - cell array containing the indices of intersecting sets
%    setType - which set has been determined to intersect the guard set
%              ('time-interval' or 'time-point')
%    Rcont - reachSet object storing the reachable set
%    options - struct with settings for reachability analysis
%
% Outputs:
%    Rguard - cell array containing the intersection with the guards
%    actGuards - cell array with indices of the active guards
%    minInd - minimum index of set intersecting guard 
%    maxInd - maximum index of set interescting guard
%
% References: 
%   [1] M. Althoff et al. "Computing Reachable Sets of Hybrid Systems Using 
%       a Combination of Zonotopes and Polytopes", 2009
%   [2] A. Girard et al. "Zonotope/Hyperplane Intersection for Hybrid 
%       Systems Reachablity Analysis"
%   [3] M. Althoff et al. "Avoiding Geometic Intersection Operations in 
%       Reachability Analysis of Hybrid Systems"
%   [4] S. Bak et al. "Time-Triggered Conversion of Guards for Reachability
%       Analysis of Hybrid Automata"
%   [5] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/reach

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       08-May-2007 
% Last update:   21-September-2007
%                30-July-2016
%                23-November-2017
%                20-April-2018 (intersect guard sets with invariant)
%                23-December-2019 (restructured the code)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check if there exist guard intersections
    relIndex = unique(setInd);

    if isempty(relIndex)
        Rguard=[]; actGuards=[]; minInd=[]; maxInd=[];
        return;
    end

    % extract the guard sets that got hit
    guardInd = unique(guards);
    Pguard = cell(length(loc.transition),1);

    for i=1:length(guardInd)
        Pguard{guardInd(i)} = loc.transition(guardInd(i)).guard;
    end

    % group the reachable sets which intersect guards
    Rtp = Rcont.timePoint.set;
    if strcmp(setType,'time-interval')
        [minInd,maxInd,P,actGuards] = aux_groupSets(Rcont.timeInterval.set,guards,setInd);
    elseif strcmp(setType,'time-point')
        [minInd,maxInd,P,actGuards] = aux_groupSets(Rtp,guards,setInd);
    end

    % loop over all guard intersections
    Rguard = cell(length(minInd),1);

    for i = 1:length(minInd)

        % get current guard set
        guard = Pguard{actGuards(i)};
        
        % remove all intersections where the flow does not point in the
        % direction of the guard set (only check if there is no instant
        % transition, i.e., there is a time-interval reachable set)
        if strcmp(setType,'time-interval') && ...
                (isa(guard,'conHyperplane') || isa(guard,'levelSet'))
            try
                [res,P{i}] = checkFlow(loc,guard,P{i},options);
                if ~res
                    continue;
                end
            end
        end

        if isa(guard,'levelSet')
            % if current guard set is a level set, only level set method
            % applies
            Rguard{i} = guardIntersect_levelSet(loc,P{i},guard);
            continue
        end

        % selected method for the calculation of the intersection
        switch options.guardIntersect
             
            % compute intersection with the method in [1]
            case 'polytope'

                Rguard{i} = guardIntersect_polytope(loc,P{i},guard,options);

            % compute intersection using constrained zonotopes
            case 'conZonotope'

                Rguard{i} = guardIntersect_conZonotope(loc,P{i},guard,options);
               
            % compute intersection with the method in [2]
            case 'zonoGirard'
                
                Rguard{i} = guardIntersect_zonoGirard(loc,P{i},guard,options);
                
            % compute intersection with the method in [3]
            case 'hyperplaneMap'
                
                R0 = aux_getInitialSet(Rtp,minInd(i));
                Rguard{i} = guardIntersect_hyperplaneMap(loc,guard,R0,options);   

            % compute intersection with the method in [4]
            case 'pancake'
                
                R0 = aux_getInitialSet(Rtp,minInd(i));
                Rguard{i} = guardIntersect_pancake(loc,R0,guard,actGuards(i),options);
                
            % compute intersection with method for nondeterministic guards
            case 'nondetGuard'
                
                Rguard{i} = guardIntersect_nondetGuard(loc,P{i},guard,options);
                
            % compute intersection with the method in [5]    
            case 'levelSet'
                
                if isa(guard,'conHyperplane')
                    guard = levelSet(guard);
                end
                Rguard{i} = guardIntersect_levelSet(loc,P{i},guard);
                
            otherwise
                throw(CORAerror('CORA:wrongFieldValue','options.guardIntersect',...
                    {'polytope','conZonotope','zonoGirard',...
                    'hyperplaneMap','pancake','nondetGuard','levelSet'}));

        end
    end
    
    % remove all empty intersections
    [Rguard,minInd,maxInd,actGuards] = aux_removeEmptySets(Rguard,minInd,...
                                                       maxInd,actGuards);
                                                   
    % convert sets back to polynomial zonotopes
    if isa(options.R0,'polyZonotope')
        for i = 1:length(Rguard)
            if isa(Rguard{i},'zonotope')
                Rguard{i} = polyZonotope(Rguard{i}); 
            end
        end
    end
    
end


% Auxiliary functions -----------------------------------------------------

function [minInd,maxInd,P,guards] = aux_groupSets(Pset,guards,setIndices)
% group the reachable sets which intersect guard sets. The sets in one
% group all intersect the same guard set and are located next to each other

    % initialization
    guardInd = unique(guards);
    setIndicesGuards = cell(length(guardInd),1);
    P = cell(length(guardInd),1);
    
    % Step 1: group according to hit guard sets 
    for i = 1:length(guardInd)
        ind = find(guards == guardInd(i));
        setIndicesGuards{i} = setIndices(ind);
        P{i} = Pset(setIndices(ind));
    end
    
    % Step 2: group accoring to the location (neigbouring sets together)
    [minInd,maxInd,P,guards] = aux_removeGaps(setIndicesGuards,guardInd,P);

end

function [minInd,maxInd,Pint,guards] = aux_removeGaps(setIndicesGuards,guards,Pint)
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

function [R,minInd,maxInd,actGuards] = aux_removeEmptySets(R,minInd,maxInd,actGuards)
% remove all sets for which the intersection with the guard set turned out
% to be empty

    % get the indices of the non-emtpy sets
    ind = ones(length(R));
    counter = 1;
    
    for i = 1:length(R)
       if isa(R{i},'zonoBundle')
           if any(cellfun(@(S) representsa_(S,'emptySet',eps),R{i}.Z))
               continue;
           end
       elseif representsa_(R{i},'emptySet',eps)
           continue;
       end
       ind(counter) = i;
       counter = counter + 1;
    end
    
    ind = ind(1:counter-1);
    
    % update variables
    R = R(ind);
    minInd = minInd(ind);
    maxInd = maxInd(ind);
    actGuards = actGuards(ind);
end

function R0 = aux_getInitialSet(Rtp,minInd)
% get the initial set

    if minInd == 1
        throw(CORAerror('CORA:specialError',...
            'The initial set already intersects the guard set!')); 
    else
        R0 = Rtp{minInd-1};
    end
end

% ------------------------------ END OF CODE ------------------------------
