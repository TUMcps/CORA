function [guards,setIndices,setType] = potInt(loc,R,options)
% potInt - determines which reachable sets potentially intersect with which
%    guard sets
%
% Syntax:
%    [guards,setIndices,setType] = potInt(loc,R,options)
%
% Inputs:
%    loc - location object
%    R - reachSet object storing reachable sets of location/reach
%    options - struct containing the algorithm settings
%
% Outputs:
%    guards - guards that are potentially intersected
%    setIndices - indices of the reachable sets that intersect the guards
%    setType - which set has been determined to intersect the guard set
%              ('time-interval' or 'time-point')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       08-May-2007 
% Last update:   26-October-2007
%                20-October-2010
%                27-July-2016
%                23-November-2017
%                03-December-2019 (NK, use approximate intersection test)
%                02-June-2023 (MW, immediate reach exit, special case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check whether time-interval solution given (if the invariant is empty,
% then the reachable set computation ends without computing any
% time-interval reachable sets and returns only the start set)
if ~isempty(R.timeInterval)
    Rset = R.timeInterval.set;
    setType = 'time-interval';
else
    Rset = R.timePoint.set;
    setType = 'time-point';
end

% number of reachable sets
nrSets = length(Rset);
% number of transitions in the location = number of guard sets
nrTrans = length(loc.transition);

% preallocate variables for output arguments (upper bound of entries)
guards = zeros(nrTrans*nrSets,1);
setIndices = zeros(nrTrans*nrSets,1);

% initialize number of intersections
counter = 1;

% loop over all guards
for i = 1:nrTrans
    
    % read out guard set and target location
    guardSet = loc.transition(i).guard;
    target = loc.transition(i).target;
    
    % only check if target location is not one of the terminal locations
    if ~all(target == options.finalLoc)
    
        % loop over all reachable sets
        for j = 1:nrSets

            % check if reachable set intersects the guard set
            if isIntersecting_(guardSet,Rset{j},'approx')
                guards(counter) = i;
                setIndices(counter) = j;
                counter = counter + 1;
            end

        end
    end
end

% remove zeros from resulting lists
guards = guards(1:counter-1);
setIndices = setIndices(1:counter-1);

% special case: when a timer is carried along as a state variable, it can
% happen that the time step size is chosen such that a time-point solution
% at t_k exactly intersects the guard set (time-triggered transition), but
% the time-interval solutions at [t_k-1,t_k], [t_k,t_k+1] are used instead

% we detect this case as follows
% - only one guard set is intersected
% - intersected guard set is a conHyperplane object
% - only one/two subsequent time-interval solutions intersect
% - time-point solution in the middle is contained in hyperplane
% - time-interval solution(s) are not contained in hyperplane
% since containment check is set-in-polytope, which uses support function
% evaluations, we only do this for sets, where the support function is
% somewhat quickly evaluated (zonotope, zonoBundle, conZonotope)

% check correct number of intersections and object classes
if ~isempty(R.timeInterval) ...
        && ( (length(guards) == 1 && length(setIndices) == 1) ...
        || ( length(guards) == 2 && guards(1) == guards(2) ...
        && length(setIndices) == 2 && diff(setIndices) == 1 ) ) ...
        && isa(loc.transition(guards(1)).guard,'conHyperplane') ...
        && ( isa(Rset{setIndices(1)},'zonotope') || isa(Rset{setIndices(1)},'zonoBundle') ...
        || isa(Rset{setIndices(1)},'conZonotope') )
    % check containment of time-point solution
    if contains_(loc.transition(guards(1)).guard,R.timePoint.set{setIndices(1)+1},'exact',1e-10)
        % ensure that time-interval solutions are not contained
        if ~any(cellfun(@(x) contains_(loc.transition(guards(1)).guard,x,'exact',1e-10),...
                Rset(setIndices),'UniformOutput',true))
            % choose time-point instead of time-interval solution(s)
            guards = guards(1);
            setIndices = setIndices(1)+1;
            setType = 'time-point';
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
