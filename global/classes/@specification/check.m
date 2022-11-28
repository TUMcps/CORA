function [res,ind] = check(spec,S,varargin)
% check - checks if a set satisfies the specification
%
% Syntax:  
%    [res,ind] = check(spec,S)
%    [res,ind] = check(spec,S,time)
%
% Inputs:
%    spec - specification object
%    S - contSet or reachSet object
%    time - time interval for the reachable set (class: interval)
%
% Outputs:
%    res - true/false whether set satisfies the specification
%    ind - index of the specification that is violated
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Author:       Niklas Kochdumper
% Written:      29-May-2020             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = true;

    % parse input arguments
    ind = []; time = [];
    
    if nargin > 2 && ~isempty(varargin{1})
        time = varargin{1}; 
    end

    % split into temporal logic and other specifications
    [spec,specLogic] = splitLogic(spec);

    if ~isempty(specLogic)
        for i = 1:size(specLogic,1)
            res = checkLogic(specLogic(i,1).set,S,time);
            if ~res
                ind = 1; return;
            end
        end
        if ~isempty(spec)
            res = check(spec,S,time);
        end
    end
    
    % check if set or single point is provided
    if isnumeric(S)
        
        % check if single or multiple points are provided
        if size(S,2) > 1
            for i = 1:size(S,2)
                [res,ind] = check(spec,S(:,i));
                if ~res
                    return; 
                end
            end
        else
            
            % loop over all specifications
            for i = 1:size(spec,1)

                % check if time frames overlap
                if isempty(time) && ~isempty(spec(i,1).time)
                    throw(CORAerror('CORA:specialError',...
                        'Timed specifications require a time interval.')); 
                end

                if isempty(spec(i,1).time) || contains(spec(i,1).time,time)

                    % different types of specifications
                    switch spec(i,1).type

                        case 'invariant'
                            res = contains(spec(i,1).set,S);

                        case 'unsafeSet'
                            res = ~contains(spec(i,1).set,S);

                        case 'safeSet'
                            res = contains(spec(i,1).set,S);

                        case 'custom'
                            res = checkCustom(spec(i,1).set,S);
                    end

                    % return as soon as one specification is violated
                    if ~res
                        ind = i; return;
                    end
                end
            end
            
        end
        
    elseif isa(S,'reachSet')

        % loop over all specifications
        for i = 1:size(spec,1)

            % loop over all reachable sets
            for k = 1:size(S,1)
                for j = 1:length(S(k).timeInterval.set)
                    [res,ind] = check(spec(i,1),S(k).timeInterval.set{j});
                    if ~res
                        ind = i; return; 
                    end
                end
            end
        end

    else

        % loop over all specifications
        for i = 1:size(spec,1)

            % check if time frames overlap
            if isempty(time) && ~isempty(spec(i,1).time)
                throw(CORAerror('CORA:specialError',...
                    'Timed specifications require a time interval.')); 
            end

            if isempty(spec(i,1).time) || isIntersecting(spec(i,1).time,time)

                % different types of specifications
                switch spec(i,1).type

                    case 'invariant'
                        res = checkInvariant(spec(i,1).set,S);

                    case 'unsafeSet'
                        res = checkUnsafeSet(spec(i,1).set,S);

                    case 'safeSet'
                        res = checkSafeSet(spec(i,1).set,S);

                    case 'custom'
                        res = checkCustom(spec(i,1).set,S);
                end

                % return as soon as one specification is violated
                if ~res
                    ind = i; return;
                end
            end
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = checkUnsafeSet(set,S)
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

function res = checkSafeSet(set,S)
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

function res = checkCustom(func,S)
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

function res = checkInvariant(set,S)
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