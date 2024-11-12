classdef verifyTime
% verifyTime - 1D interval representing the times at which a specification
%    must be checked (verified/falsified); crucially, the interval may
%    consist of multiple disjoint partial intervals
%    note: originally implemented for the verification algorithm [1, Alg. 3]
%    to speed up computations (using multiple intervals causes overhead)
%
% Syntax:
%    Itime = verifyTime(tStart,tFinal)
%
% Inputs:
%    tStart - start of time horizon
%    tFinal - end of time horizon
%    
% Outputs:
%    Itime - generated verifyTime object
%
% Example:
%    Itime = verifyTime([0,2; 4,5]);
%
% References:
%    [1] M. Wetzlinger et al. "Fully automated verification of linear
%        systems using inner-and outer-approximations of reachable sets",
%        TAC, 2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/private/reach_adaptive

% Authors:       Mark Wetzlinger
% Written:       10-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = public, GetAccess = public)  % Access = private
    
    bounds;     % sx2 array: sequence of s pair-wise lower and upper bounds
                % example with three disjoint parts: [0,1; 1.5,2; 3,4]

end

methods
    % constructor
    function Itime = verifyTime(varargin)

        narginchk(0,1);
        % 0. empty call
        if nargin == 0
            Itime.bounds = zeros(0,2);
            return
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'verifyTime')
            Itime = verifyTime(varargin{1}.bounds); return
        end

        % 2. parse input arguments: varargin -> vars
        bounds = varargin{1};

        % 3. check correctness of input arguments
        aux_checkInputArgs(bounds);

        % 4. assign properties
        Itime.bounds = bounds;
    end

    function Itimes = shift(Itimes,t)
        % shift time intervals by t (e.g., if tStart not 0)
        for i=1:numel(Itimes)
            Itimes(i).bounds = Itimes(i).bounds + t;
        end
    end

    function tStart = startTime(Itimes)
        % read out minimum start time over all intervals
        allBounds = vertcat(Itimes.bounds);
        tStart = min(allBounds(:,1));
    end

    function tFinal = finalTime(Itimes)
        % read out maximum final time over all intervals
        allBounds = vertcat(Itimes.bounds);
        tFinal = max(allBounds(:,2));
    end

    function s = numIntervals(Itimes)
        % read out number of disjoint time intervals
        s = arrayfun(@(Itime) size(Itime.bounds,1),...
                     Itimes,'UniformOutput',true);
    end

    function I = interval(Itime,varargin)
        % convert to interval
        narginchk(1,2);

        % by default: exact conversion (only possible for single interval!)
        mode = setDefaultValues({'exact'},varargin);
        inputArgsCheck({{Itime,'att','verifyTime','scalar'};...
                        {mode,'str',{'exact','outer'}}});

        switch mode
            case 'exact'
                Itime = compact(Itime);
                if numIntervals(Itime) > 1
                    throw(CORAerror('CORA:noExactAlg',...
                        ['Multiple disjoint time intervals cannot be '...
                        'exactly converted into a single interval.']));
                end
                I = interval(Itime.startTime(),Itime.finalTime());
            case 'outer'
                I = interval(Itime.startTime(),Itime.finalTime());
        end

    end

    function Itime = compact(Itime,varargin)
        % unify adjacent intervals (if possible)
        narginchk(1,2);

        % quick exit: already in minimal form (empty, single time interval)
        s = numIntervals(Itime);
        if s <= 1
            return
        end

        % set tolerance
        tol = setDefaultValues({eps*(Itime.finalTime() - Itime.startTime())},varargin);

        % indices where the end matches the next start
        idx = withinTol(Itime.bounds(1:end-1,2), Itime.bounds(2:end,1),tol);

        % quick exits: all match/no matches
        if all(idx)
            Itime.bounds = [Itime.startTime(), Itime.finalTime()];
            return
        elseif ~any(idx)
            return
        end
        
        % loop over all intervals
        idx = [false;idx]; i = 1;
        while i <= s
            % for all rows that should be removed, see if there are
            % multiple rows in succession that should be removed
            if idx(i)
                untilIdx = find(~idx(i+1:end),1,'first') - 1;
                if isempty(untilIdx)
                    untilIdx = 0;
                end
                lastIdx = i + untilIdx;
                Itime.bounds(i-1,2) = Itime.bounds(lastIdx,2);
                i = lastIdx;
            end
            i = i + 1;
        end
        % remove rows which have been integrated into upper rows
        Itime.bounds(idx,:) = [];

    end

    function res = representsa(Itimes,type,tol)
        % check if verifyTime object represents a contSet class
        if strcmp(type,'emptySet')
            res = arrayfun(@(Itime) size(Itime.bounds,1) == 0,...
                           Itimes,'UniformOutput',true);
        elseif strcmp(type,'interval')
            Itimes = compact(Itimes,tol);
            res = arrayfun(@(Itime) size(Itime.bounds,1) == 1,...
                           Itimes,'UniformOutput',true);
        else
            throw(CORAerror('CORA:notSupported',...
                'Only comparison to emptySet and interval supported.'));
        end
    end

    function res = isequal(Itime,S)
        % comparison only implemented for verifyTime objects
        inputArgsCheck({{Itime,'att','verifyTime'};...
                        {S,'att','verifyTime'}});

        % case where one or both are empty
        Itime_empty = representsa(Itime,'emptySet',0);
        S_empty = representsa(S,'emptySet',0);
        if xor(Itime_empty,S_empty)
            res = false;
            return
        elseif Itime_empty && S_empty
            res = true;
            return
        end

        % rewrite both in minimal representation
        Itime = compact(Itime);
        S = compact(S);

        % compare number of time intervals
        if numIntervals(Itime) ~= numIntervals(S)
            res = false; 
            return
        end

        % compare bounds up to tolerance (sizes must match since we checked
        % the number of disjoint time intervals above)
        tol = eps*(Itime.finalTime() - Itime.startTime());
        res = compareMatrices(Itime.bounds,S.bounds,tol,"equal",true);
    end
    function res = eq(Itime,S)
        res = isequal(Itime,S);
    end
    function res = ne(Itime,S)
        res = ~isequal(Itime,S);
    end

    function [res,idx] = contains(Itime,t)
        % determine whether a disjoint time interval contains a specific
        % point in time; we also return the index of the time interval
        % where the point in time is contained (otherwise [])
        if representsa(Itime,'emptySet',0)
            res = false; idx = [];
            return
        end

        % must be larger than lower bound and smaller than upper bound
        logIdx = Itime.bounds(:,1) <= t & Itime.bounds(:,2) >= t;
        res = any(logIdx);
        idx = find(logIdx,1,'first');
    end

    function res = isIntersecting(Itime,tInt)
        % check if a time interval intersects any (part) of the disjoint
        % intervals
        if representsa(Itime,'emptySet',0)
            res = false;
            return
        end

        % check for intersection
        res = ~( tInt(2) <= Itime.startTime() || tInt(1) > Itime.finalTime() || ...
            any(tInt(1) >= Itime.bounds(1:end-1,2) & tInt(2) <= Itime.bounds(2:end,1)));
    end

    function [tSwitch,inside] = timeUntilSwitch(Itime,t)
        % determines the time from a given point in time until the end of 
        % the current time interval or the start of the next interval;
        % additionally, we return whether we continue inside (true) the
        % intervals or outside (false)
        if representsa(Itime,'emptySet',0)
            tSwitch = []; inside = [];
            return
        end

        % compact representation
        Itime = compact(Itime);

        allTimes = reshape(Itime.bounds',[],1);
        idx = find(t < allTimes & ~withinTol(t,allTimes),1,'first');
        tSwitch = allTimes(idx) - t;
        % switch goes from outside to inside if 
        inside = mod(idx,2) == 0;
    end

    function Itime = setdiff(Itime,tInt)
        % remove time interval from verifyTime object

        Itime.bounds = Itime.bounds;
        set_unsat_col = reshape(Itime.bounds',numel(Itime.bounds),1);
        if mod(sum(tInt(1) >= set_unsat_col),2) ~= 0
            % lower bound starts inside unverified time interval
            % t(1) \in [ timeInterval )
            idx = find(tInt(1) >= Itime.bounds(:,1) & tInt(1) <= Itime.bounds(:,2));
            
            if tInt(2) <= Itime.bounds(idx,2)
                if tInt(1) > Itime.bounds(idx,1)
                    Itime.bounds = [Itime.bounds(1:idx-1,:); [Itime.bounds(idx,1), tInt(1)]; Itime.bounds(idx:end,:)];
                    idx = idx + 1;
                end
                % split, potential merge later
                Itime.bounds(idx,1) = tInt(2);
                tInt = [];
            else
                % remove interval, potential merge later
                Itime.bounds(idx,2) = tInt(1);
                if idx < size(Itime.bounds,1)
                    tInt(1) = Itime.bounds(idx+1,1);
                end
                if tInt(2) <= tInt(1)
                    tInt = [];
                end
            end
        end
        
        % now: lower bound starts in between unverified time intervals or at
        % maximum at the start point of an unverified set
        % t(1) \in [ notTimeInterval )
        while ~isempty(tInt)
            idx = find(tInt(1) <= Itime.bounds(:,1),1,'first');
            % upper bound is at least at the beginning of the next time interval
            if tInt(2) <= Itime.bounds(idx,2)
                % split, potential merge later
                Itime.bounds(idx,1) = tInt(2);
                tInt = [];
            else
                % remove entire thing (full time interval verified)
                if idx < size(Itime.bounds,1)
                    tInt(1) = Itime.bounds(idx,2);
                    if tInt(2) < Itime.bounds(idx+1,1)
                        tInt = [];
                    end
                else
                    tInt = [];
                end
                Itime.bounds(idx,:) = [];
            end
        end
        
        % remove
        idxRemove = abs(Itime.bounds(:,2) - Itime.bounds(:,1)) < 1e-14;
        Itime.bounds(idxRemove,:) = [];

    end

    function Itime = unify(Itimes)
        % unifies the intervals of multiple verifyTime objects into a
        % single verifyTime object

        % initialize bounds of resulting object
        Itime = Itimes(1);
        % read out bounds
        bnds = Itime.bounds;

        % read out start time and final time over all objects, convert to
        % interval for quick exit check below
        Itime_full = verifyTime([min(Itimes.startTime()),...
                                 max(Itimes.finalTime())]);

        % loop over all other objects
        for i=2:numel(Itimes)
        
            % quick exit if maximum possible time interval already covered
            if numIntervals(Itime) == 1 && isequal(Itime,Itime_full)
                break
            end
            
            % loop over all time intervals
            for j=1:numIntervals(Itimes(i))
                
                % expand the j-th time interval for better readability
                t = Itimes(i).bounds(j,:);
                
                % fully inside
                idx = find(t(1) >= bnds(:,1) & t(2) <= bnds(:,2),1);
                if ~isempty(idx)
                    % already entirely covered -> continue
                    continue
                end

                % cases with overlap: cut overlap
                idx = find(t(1) >= bnds(:,1) & t(1) <= bnds(:,2),1);
                if ~isempty(idx)
                    t(1) = bnds(idx,2);
                end
                idx = find(t(2) >= bnds(:,1) & t(2) <= bnds(:,2),1);
                if ~isempty(idx)
                    t(2) = bnds(idx,1);
                end

                % check if current time interval contains an existing time interval
                idx = find(t(1) <= bnds(:,1) & t(2) >= bnds(:,2));
                if ~isempty(idx)
                    bnds(idx,:) = [];
                    % will be implicitly re-added later as t contains it
                    if isempty(bnds)
                        bnds = t;
                        continue
                    end
                end
        
        
                % simple cases: no overlap (maximum: point-wise)
                if t(2) < bnds(1,1)
                    bnds = [t; bnds];
                elseif t(2) == bnds(1,1)
                    bnds(1,1) = t(1);
                elseif t(1) == bnds(end,2)
                    bnds(end,2) = t(2);
                elseif t(1) > bnds(end,2)
                    bnds = [bnds; t];
                else
                    idx = find(t(1) >= bnds(1:end-1,2) & t(2) <= bnds(2:end,1));
                    if ~isempty(idx)
                        bnds = [bnds(1:idx,:); t; bnds(idx+1:end,:)];
                        idxMerge = 0;
                        while ~isempty(idxMerge)
                            idxMerge = find(bnds(2:end,1) == bnds(1:end-1,2));
                            if ~isempty(idxMerge)
                                bnds = [bnds(1:idx-1,:);
                                             [bnds(idx,1), bnds(idx+1,2)];
                                             bnds(idx+2:end,:)];
                            end
                        end
                    end
                end
                
            end            
        end

        % assign new bounds to resulting object
        Itime.bounds = bnds;
    end

end

end


% Auxiliary functions -----------------------------------------------------

function aux_checkInputArgs(bounds)
    % check input arguments

    % size must be sx2
    if ~isnumeric(bounds) || size(bounds,2) ~= 2
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Bounds must be a s-by-2 numeric array'));
    end

    % values must be finite
    if ~all(isfinite(bounds),'all')
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Bounds must contain finite values (no Inf or NaN)'));
    end

    % bounds must be proper
    if any(bounds(:,1) > bounds(:,2))
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Time intervals must be proper (lower bounds <= upper bounds'));
    end

    % time intervals must be in order (start of next interval must be later
    % than end of previous interval)
    if any(bounds(1:end-1,2) > bounds(2:end,1))
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Time intervals must be proper (lower bounds <= upper bounds'));
    end

end

% ------------------------------ END OF CODE ------------------------------
