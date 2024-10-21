function assertLoop(cond,varargin)
% assertLoop - this function checks if <cond> is satisfied within a loop
%    and adds the given seed and for loop counters to the assert message.
%    Note: If the counter variables cannot be determined, they are shown as 'c<i>=%i'.
%
% Syntax:
%    res = assertLoop(cond,varargin)
%    res = assertLoop(cond,msg,seed,varargin)
%
% Inputs:
%    cond - logical
%    msg - additional message
%    seed - seed of current iteration
%    varargin - loop counters
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: assert

% Authors:       Tobias Ladner
% Written:       29-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no early exit if cond=true as it is important to check whether something
% fails in case cond=false.

% parse input
narginchk(1,Inf);
[cond,msg,seed,loopCounters,counterOffset] = aux_parseInput(cond,varargin);

% build message
assertmsg = 'Assertion failed.';

if ~isempty(msg)
    assertmsg = sprintf('%s %s',assertmsg,msg);
    if ~assertmsg(end) == '.'
        assertmsg = sprintf('%s.',assertmsg);
    end
end

% add loop counters
if ~isempty(loopCounters)
    loopCountersText = cell(1,numel(loopCounters));
    for i = 1:numel(loopCounters)
        counterName = inputname(counterOffset+i);
        if isempty(counterName)
            % name is empty even if any previous inputs are weird...
            counterName = sprintf('c%i',i);
        end
        if isnumeric(loopCounters{i})
            loopCountersText{i} = sprintf('%s=%i',counterName,loopCounters{i});
        else
            loopCountersText{i} = sprintf('%s=''%s''',counterName,loopCounters{i});
        end
    end
    % joining with 'j' for easy copy-paste into command window
    countersmsg = strjoin(loopCountersText,'; ');
    if isscalar(loopCountersText)
        assertmsg = sprintf('%s Loop counter: %s;', assertmsg, countersmsg);
    else
        assertmsg = sprintf('%s Loop counters: %s;', assertmsg, countersmsg);
    end
end

% add seed
if ~isempty(seed)
    assertmsg = sprintf('%s seed=%i;',assertmsg,seed);
end

% do assertion
assert(all(cond,'all'),assertmsg);

end


% Auxiliary functions -----------------------------------------------------

function [cond,msg,seed,loopCounters,counterOffset] = aux_parseInput(cond,givenvalues)

    % parse input
    [msg,seed] = setDefaultValues({'',[]},givenvalues);
    if isnumeric(msg)
        % no msg and seed given, all given values are loop counters
        msg = ''; seed = [];
        loopCounters = givenvalues;
        counterOffset = 1;
    else
        loopCounters = givenvalues(3:end);
        counterOffset = 3;
    end

    % check input
    inputArgsCheck({ ...
        {cond,'att','logical'}, ...
        {msg,'att',{'char','string'}}, ...
        {seed,'att','numeric'}, ...
        {loopCounters,'att','cell'}, ...
    });

    % additional checks
    if ~isscalar(seed) && ~isempty(seed)
        throw(CORAerror('CORA:wrongValue','third','Seed has to be a numeric scalar or empty.'));
    end
    if isempty(loopCounters)
       throw(CORAerror('CORA:wrongValue','variable','At least one loop counter should be specified.'));
    end
    if ~all(cellfun(@(loopCounter) (isnumeric(loopCounter) && isscalar(loopCounter)) || ischar(loopCounter) || isstring(loopCounter),loopCounters))
       throw(CORAerror('CORA:wrongValue','variable','Loop counters have to be scalar integers.'));
    end
end

% ------------------------------ END OF CODE ------------------------------
