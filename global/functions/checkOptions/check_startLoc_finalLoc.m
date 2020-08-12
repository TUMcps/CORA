function options = check_startLoc_finalLoc(options, obj, ispha)
% check_startLoc_finalLoc - checks if options.startLoc/finalLoc
%  1) exist
%  2) take an allowed value
%
% Syntax:
%    options = check_startLoc_finalLoc(options, obj, 1)
%
% Inputs:
%    options - options for object
%    obj     - hybrid automaton obj (flat or parallel)
%    ispha   - flag if hybrid automaton is flat (0) or parallel (1)
%
% Outputs:
%    options   - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      05-Mar-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

% handle ispha with nargin...

strct = 'params';
option1 = 'startLoc';
option2 = 'finalLoc';
if ispha == 0
    % flat hybrid automaton
    locations = obj.location;
    numLoc = length(locations);
    
    % startLoc has to be int and < numLoc (without 0)
    if ~isfield(options,option1)
        error(printOptionMissing(obj,option1,strct))
    elseif options.startLoc > numLoc || mod(options.startLoc,1.0) ~= 0
        if options.startLoc <= 0
            error(printOptionOutOfRange(obj,option1,strct));
        end
    end

    % finalLoc has to be int and < numLoc
    if ~isfield(options,option2)
        options.finalLoc = 0; % no final location as default
    elseif options.finalLoc > numLoc || mod(options.finalLoc,1.0) ~= 0
        if options.finalLoc < 0
            error(printOptionOutOfRange(obj,option2,strct));
        end
    end
    
elseif ispha == 1
    % parallel hybrid automaton 
    numComp = length(obj.components);
    
    % startLoc has to be a vector with as much rows as number of components
    % feasible values dependent on number of locations in each component
    if ~isfield(options,option1)
        error(printOptionMissing(obj,option1,strct));
    elseif ~all(size(options.startLoc) == [numComp,1])
        error(printOptionOutOfRange(obj,option1,strct)); 
    else
        for c=1:numComp
            % loop through every component
            comp = obj.components{c};
            numLoc = length(comp.location);
            if options.startLoc(c) > numLoc || options.startLoc(c) <= 0 || ...
               mod(options.startLoc(c),1.0) ~= 0 
                    error(printOptionOutOfRange(obj,option1,strct));
            end
        end
    end

    % finalLoc has to be a vector with as much rows as number of components
    % feasible values dependent on number of locations in each component
    if ~isfield(options,option2)
        options.finalLoc = zeros(numComp,1);
    elseif ~all(size(options.finalLoc) == [numComp,1])
        error(printOptionOutOfRange(obj,option2,strct)); 
    else
        for c=1:numComp
            % loop through every component
            comp = obj.components{c};
            numLoc = length(comp.location);
            if options.finalLoc(c) > numLoc || options.finalLoc(c) < 0 || ...
               mod(options.finalLoc(c),1.0) ~= 0 
                    error(printOptionOutOfRange(obj,option2,strct));
            end
        end
    end

end


end

%------------- END OF CODE --------------
