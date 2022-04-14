function [Rguard,actGuards,minInd,maxInd] = guardIntersect(obj,guards,setInd,Rcont,options)
% guardIntersect - computes and enclosure of the intersection between the
%                  reachable set and the guard sets
%
% Syntax:  
%    [Rguard,actGuards,minInd,maxInd] = 
%                           guardIntersect(obj,guards,setInd,Rcont,options)
%
% Inputs:
%    obj - location object
%    guards - cell array containing the guard sets that have been hit
%    setInd - cell array containing the indices of intersecting sets
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      08-May-2007 
% Last update:  21-September-2007
%               30-July-2016
%               23-November-2017
%               20-April-2018 (intersect guard sets with invariant)
%               23-December-2019 (restructured the code)
% Last revision:---

%------------- BEGIN CODE --------------

    % check if the there exist guard intersections
    relIndex = unique(setInd);

    if isempty(relIndex)
        Rguard=[]; actGuards=[]; minInd=[]; maxInd=[];
        return;
    end

    % extract the guard sets that got hit
    guardInd = unique(guards);
    Pguard = cell(length(obj.transition),1);

    for i=1:length(guardInd)
        Pguard{guardInd(i)} = obj.transition{guardInd(i)}.guard;
    end

    % extract time interval and time point reachable set
    R = Rcont.timeInterval.set;
    Rtp = Rcont.timePoint.set;

    % group the reachable sets which intersect guards
    [minInd,maxInd,P,actGuards] = groupSets(R,guards,setInd);
    

    % loop over all guard intersections
    Rguard = cell(length(minInd),1);

    for i = 1:length(minInd)

        % get current guard set
        guard = Pguard{actGuards(i)};
        
        % remove all intersections where the flow does not point in the
        % direction of the guard set
        if isa(guard,'conHyperplane') || isa(guard,'levelSet')
            [res,P{i}] = checkFlow(obj,guard,P{i},options);
            if ~res
                continue;
            end
        end

        % selected method for the calculation of the intersection
        switch options.guardIntersect
             
            % compute intersection with the method in [1]
            case 'polytope'

                Rguard{i} = guardIntersect_polytope(obj,P{i},guard,options);

            % compute intersection using constrained zonotopes
            case 'conZonotope'

                Rguard{i} = guardIntersect_conZonotope(obj,P{i},guard,options);
               
            % compute intersection with the method in [2]
            case 'zonoGirard'
                
                Rguard{i} = guardIntersect_zonoGirard(obj,P{i},guard,options);
                
            % compute intersection with the method in [3]
            case 'hyperplaneMap'
                
                R0 = getInitialSet(Rtp,minInd(i));
                Rguard{i} = guardIntersect_hyperplaneMap(obj,guard,R0,options);   

            % compute intersection with the method in [4]
            case 'pancake'
                
                R0 = getInitialSet(Rtp,minInd(i));
                Rguard{i} = guardIntersect_pancake(obj,R0,guard,actGuards(i),options);
                
            % compute intersection with method for nondeterministic guards
            case 'nondetGuard'
                
                Rguard{i} = guardIntersect_nondetGuard(obj,P{i},guard,options);
                
            % compute intersection with the metohd in [5]    
            case 'levelSet'
                
                Rguard{i} = guardIntersect_levelSet(obj,P{i},guard);
                
            otherwise
                error('Wrong value for "options.guardIntersect"!');

        end
    end
    
    % remove all empty intersections
    [Rguard,minInd,maxInd,actGuards] = removeEmptySets(Rguard,minInd, ...
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



% Auxiliary Functions -----------------------------------------------------

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

function [R,minInd,maxInd,actGuards] = removeEmptySets(R,minInd,maxInd,actGuards)
% remove all sets for which the intersection with the guard set turned out
% to be empty

    % get the indices of the non-emtpy sets
    ind = ones(length(R));
    counter = 1;
    
    for i = 1:length(R)
       if isa(R{i},'zonoBundle')
           if any(cellfun(@isempty,R{i}.Z))
               continue;
           end
       elseif isempty(R{i})
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

function R0 = getInitialSet(Rtp,minInd)
% get the initial set

    if minInd == 1
       error('The initial set already intersects the guard set!'); 
    else
       R0 = Rtp{minInd-1};
    end
end

%------------- END OF CODE --------------