function sig = evaluateSignal(phi, dur, signals, logic)
% evaluateSignal - evalutate the STL formula on the given signals
%
% Note that this method only supports desugared formulas.
%
% Syntax:
%    sig = evaluateSignal(phi, dur, signals);
%
% Inputs:
%    phi - STL formula
%    dur - duration of the resulting signal
%    signals - map of atomic propositions to their signal
%    logic - which logic to use (kleene or boolean)
%
% Outputs:
%    sig - signal describing the validity of the formula
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl, signal

% Authors:       Benedikt Seidl
% Written:       24-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    arguments
        phi stl
        dur double {mustBeNonnegative}
        signals containers.Map
        logic = 'kleene'
    end
    
    inputDuration = dur + maximumTime(phi);
    
    sigs = values(signals);
    
    for i = 1:length(sigs)
        % check length of all signals
        if duration(sigs{i}) < inputDuration
            throw(CORAerror('CORA:wrongValue', 'signals/dur', ...
                ['The signals must be as long as the duration of the ' ...
                 'requested signal plus the horizon of the formula.']));
        end
    
        % check type of all signals
        if ~strcmp(class(sigs{i}.value), logic)
            throw(CORAerror('CORA:wrongValue', 'signals/logic', ...
                'The signals must have the correct logic.'));
        end
    end
    
    % check signals of atomic propositions
    % if ~all(isKey(signals, values(propositions(phi))))
    %     throw(CORAerror('CORA:wrongValue', 'signals/phi', ...
    %         'All atomic propositions must be defined'));
    % end
    
    % set default logical values
    if strcmp(logic, 'logical')
        tt = true;
        ff = false;
    elseif strcmp(logic, 'kleene')
        tt = kleene.True;
        ff = kleene.False;
    end
    
    % call helper function
    sig = cutoff(inner(phi), dur);
    
    function sig = inner(phi)
        % pattern matching on the formula
        if strcmp(phi.type, 'true')
            sig = signal(inputDuration, tt);
        elseif strcmp(phi.type, 'false')
            sig = signal(inputDuration, ff);
        elseif strcmp(phi.type, 'variable')
            sig = cutoff(signals(formattedDisplayText(phi)), inputDuration);
        elseif strcmp(phi.type, '&')
            sig = inner(phi.lhs) & inner(phi.rhs);
        elseif strcmp(phi.type, '|')
            sig = inner(phi.lhs) | inner(phi.rhs);
        elseif strcmp(phi.type, '~')
            sig = ~ inner(phi.lhs);
        elseif strcmp(phi.type, 'until')
            lsig = inner(phi.lhs);
            rsig = inner(phi.rhs);
    
            int = interval(phi.from, phi.to);
    
            if strcmp(logic, 'logical')
                % calculate until operator on signal
                sig = until(lsig, int, rsig);
            elseif strcmp(logic, 'kleene')
                % calculate until operator for unknown and true regions
                sigU = until(lsig, int, rsig, kleene.False, kleene.Unknown);
                sigT = until(lsig, int, rsig, kleene.False, kleene.True);
    
                % combine both regions
                sig = sigU | sigT;
            end
        else
            throw(CORAerror('CORA:notSupported',...
                  'Only desugared formulas are supported.'));
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
