function res = checkLivelock(tracker)
% checkLivelock - checks for livelock in current analysis of parallel
%    hybrid automaton
% note: currently only designed for one branch in the reachable set
%
% Syntax:
%    res = checkLivelock(tracker)
%
% Inputs:
%    tracker (struct) - meta data about the current run, with fields
%       .switchingTime: time when locations change
%       .locID: location identifiers for composed automaton
%       .transitions: transitions that each component has taken
%           (0 denotes a virtual self-transition to facilitate the merge of
%           transitions -- not part of the original list of transitions)
%       .syncLabel: synchronization label assigned to each transition
%
% Outputs:
%    res - true/false whether livelock has occurred
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton/reach

% Authors:       Mark Wetzlinger
% Written:       01-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume no livelock
res = false;

% number of iterations in main loop until now
numIter = length(tracker);

% quick exit: first entry into function
if numIter == 1
    return
end

% quick exit: time has passed since between last two locIDs
if ~isequal(tracker(end-1).switchingTime,tracker(end).switchingTime)
    return
end

% procedure: starting from last entry, list all locIDs that have occurred
% at the same time; if there is a repetition -- including repeated composed
% transitions and their associated synchronization labels -- then a
% livelock has been detected

% count backwards
for i=numIter-1:-1:1
    if isequal(tracker(i).switchingTime,tracker(end).switchingTime) ...
            
        if all(tracker(i).locID == tracker(end).locID) ...
            && all(tracker(i).transition == tracker(end).transition) ...
            && strcmp(tracker(i).syncLabel,tracker(end).syncLabel)

            % livelock occurred: print message on command window
            fprintf("\n");
            disp("Livelock detected at time t=" + ...
                string(tracker(end).switchingTime) + ".");
            disp("The reachable sets until the current step are returned.");
            fprintf("\n");
    
            % set flag to true so that main loop can exit
            res = true;
            break

        end

    else
        % stop as no other indices can be true (time only moves forward)
        break
    end
end

% ------------------------------ END OF CODE ------------------------------
