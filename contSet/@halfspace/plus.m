function hs = plus(summand1,summand2)
% plus - Overloaded '+' operator for the addition of a vector with a
%    halfspace
%
% Syntax:  
%    hs = plus(summand1,summand2)
%
% Inputs:
%    summand1 - halfspace object or numerical vector
%    summand2 - halfspace object or numerical vector
%
% Outputs:
%    hs - halfspace object
%
% Example: 
%    hs = halfspace([1 1],2);
%    summand = [1; -0.5];
%    hs + summand
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      28-August-2013
% Last update:  16-March-2021 (MW, error handling)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing: assign halfspace and summand
[hs,summand] = findClassArg(summand1,summand2,'halfspace');

try

    % compute Minkowski sum
    if isnumeric(summand)
        hs.d = hs.d + hs.c.'*summand;
    else
        % no other summands are currently implemented
        throw(CORAerror('CORA:noops',hs,summand));
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if isempty(hs)
        return
    elseif (isnumeric(summand) && isempty(summand)) ...
            || (isa(summand,'contSet') && isemptyobject(summand))
        hs = halfspace(); return
    end

    % check whether different dimension of ambient space
    equalDimCheck(hs,summand);

    % other error...
    rethrow(ME);

end

%------------- END OF CODE --------------