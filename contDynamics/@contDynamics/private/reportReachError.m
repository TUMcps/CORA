function reportReachError(ME,time,ind)
% reportReachError - reports information to the user in case
%    the reachable set explodes; in other cases, the same information
%    as without this function is displayed
%
% Syntax:  
%    reportReachError(ME,time,ind)
%
% Inputs:
%    ME - MException object
%    time - current time
%    ind - cuurnend step
%
% Outputs:
%    ---

% Author:       Mark Wetzlinger
% Written:      19-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % called due to set explosion
    if strcmp(ME.identifier,'reach:setexplosion')

        % display information to user
        fprintf("\n");
        disp(ME.message);
        disp("  Step " + ind + " at time t=" + time);
        disp("The reachable sets until the current step are returned.");
        fprintf("\n");
    else
        % any other run-time error: report information
        rethrow(ME);
    end
end

%------------- END OF CODE --------------