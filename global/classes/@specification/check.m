function [res,indSpec,indObj] = check(spec,S,varargin)
% check - checks if a set satisfies the specification
%
% Syntax:  
%    [res,indSpec,indObj] = check(spec,S)
%    [res,indSpec,indObj] = check(spec,S,time)
%
% Inputs:
%    spec - specification object
%    S - numeric, contSet, reachSet, or simResult object
%    time - time interval for the reachable set (class: interval)
%
% Outputs:
%    res - true/false whether set satisfies the specification
%    indSpec - index of the first specification that is violated
%    indObj - index of the first object S that is violated
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      29-May-2020             
% Last update:  22-March-2022 (TL: simResult, indObj)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin < 2
    throw(CORAerror("CORA:notEnoughInputArgs", 2))
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs', 3))
end
time = setDefaultValues({[]}, varargin);

% init outputs
res = true; indSpec = []; indObj = [];

% start checking ----------------------------------------------------------

% split into temporal logic and other specifications
[spec,specLogic] = splitLogic(spec);

% check temporal logic ----------------------------------------------------

if ~isempty(specLogic)
    for i = 1:size(specLogic,1)
        res = checkLogic(specLogic(i,1).set,S,time);
        if ~res
            indSpec = 1; return;
        end
    end
    if ~isempty(spec)
        res = check(spec,S,time);
    end
end

% check other specification -----------------------------------------------

% check S object
if isnumeric(S) % ---------------------------------------------------------
    
    % check if single or multiple points are provided
    if size(S,2) > 1
        % multiple points
        for i = 1:size(S,2)
            [res,indSpec] = check(spec,S(:,i),time(i));
            if ~res
                indObj = i;
                return; 
            end
        end

    else
        % single point
        
        % loop over all specifications
        for i = 1:size(spec,1)
            spec_i = spec(i);

            % check if time frames overlap
            if isempty(time) && ~isempty(spec_i.time)
                throw(CORAerror('CORA:specialError',...
                    'Timed specifications require a time interval.')); 
            end

            if isempty(spec_i.time) || contains(spec_i.time,time)

                % different types of specifications
                switch spec_i.type

                    case 'invariant'
                        res = contains(spec_i.set,S);

                    case 'unsafeSet'
                        res = ~contains(spec_i.set,S);

                    case 'safeSet'
                        res = contains(spec_i.set,S);

                    case 'custom'
                        res = aux_checkCustom(spec_i.set,S);
                end

                % return as soon as one specification is violated
                if ~res
                    indSpec = i; 
                    indObj = 1;
                    return;
                end
            end
        end
    end
    
elseif isa(S,'simResult') % -----------------------------------------------

    % loop over all simulations
    for i = 1:length(S)
        S_i = S(i);
        for k = 1:length(S_i.x)
            [res,indSpec] = check(spec, S_i.x{k}', interval(S_i.t{k}));
            if ~res
                indObj = {i,k};
                return; 
            end
        end
    end
    
elseif isa(S,'reachSet') % ------------------------------------------------

    % loop over all specifications
    for i = 1:size(spec,1)
        spec_i = spec(i);

        % loop over all reachable sets
        for k = 1:size(S,1)

            timePoint = S(k).timePoint;
            timeInterval = S(k).timeInterval;

            for j = 1:length(timeInterval.set)
                % where is the specification active?

                if isempty(spec_i.time) % entire time horizon
                    res = check(spec_i, timeInterval.set{j});

                elseif rad(spec_i.time) > 0 % part of time horizon
                    res = check(spec_i, ...
                        timeInterval.set{j}, ...
                        timeInterval.time{j});

                else % only active in one time point
                    % check if its on boundary of timeInterval
                    if contains(spec_i.time, timePoint.time{j})
                        % only start of time interval
                        res = check(spec_i, ...
                            timePoint.set{j}, ...
                            interval(timePoint.time{j}));

                    elseif contains(spec_i.time, timePoint.time{j+1})
                        % only end of time interval
                        res = check(spec_i, ...
                            timePoint.set{j+1}, ...
                            interval(timePoint.time{j+1}));

                    else % intermediate
                        res = check(spec_i, ...
                            timeInterval.set{j}, ...
                            timeInterval.time{j});
                    end
                end

                if ~res
                    indSpec = i; 
                    indObj = {k,j};
                    return; 
                end
            end
        end
    end

else % contSet ------------------------------------------------------------

    % loop over all specifications
    for i = 1:size(spec,1)
        spec_i = spec(i);

        % check if time frames overlap
        if isempty(time) && ~isempty(spec_i.time)
            throw(CORAerror('CORA:specialError',...
                'Timed specifications require a time interval.')); 
        end

        if isempty(spec_i.time) || isIntersecting(spec_i.time,time)

            % different types of specifications
            switch spec_i.type

                case 'invariant'
                    res = aux_checkInvariant(spec_i.set,S);

                case 'unsafeSet'
                    res = aux_checkUnsafeSet(spec_i.set,S);

                case 'safeSet'
                    res = aux_checkSafeSet(spec_i.set,S);

                case 'custom'
                    res = aux_checkCustom(spec_i.set,S);
            end

            % return as soon as one specification is violated
            if ~res
                indSpec = i; 
                indObj = 1;
                return;
            end
        end
    end
end

end

% Auxiliary Functions -----------------------------------------------------

function res = aux_checkUnsafeSet(set,S)
% check if reachable set intersects the unsafe sets

    if iscell(S)
        res = true;
        for i = 1:length(S)
            try
                res = ~isIntersecting(set,S{i}); 
            catch
                res = ~isIntersecting(set,S{i},'approx'); 
            end
            if ~res
               return; 
            end
        end   
    else
        try
            res = ~isIntersecting(set,S);
        catch
            res = ~isIntersecting(set,S,'approx'); 
        end
    end
end

function res = aux_checkSafeSet(set,S)
% check if reachable set is inside the safe set

    if iscell(S)
        res = true;
        for i = 1:length(S)
           res = contains(set,S{i}); 
           if ~res
              return; 
           end
        end   
    else
        res = contains(set,S);
    end
end

function res = aux_checkCustom(func,S)
% check if the reachable set satisfies a user provided specification

    if iscell(S)
        res = false;
        for i = 1:length(S)
            res = func(S{i});
            if res
               return; 
            end
        end
    else
        res = func(S);
    end
end

function res = aux_checkInvariant(set,S)
% check if reachable set intersects the invariant

    if iscell(S)
        res = false;
        for i = 1:length(S)
            res = isIntersecting(set,S{i},'approx');
            if res
               return; 
            end
        end
    else
        res = isIntersecting(set,S,'approx');
    end
end

function res = checkLogic(set,S,time)
% check if the reachable set satisfies a temporal logic specification

    if isnumeric(S)
        res = modelCheckTrace(set,S,time);
    else
        res = modelCheckReachSet(set,S);
    end
end

%------------- END OF CODE --------------