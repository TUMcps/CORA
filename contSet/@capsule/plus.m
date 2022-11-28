function C = plus(summand1,summand2)
% plus - Overloaded '+' operator for the over-approximative Minkowski 
%    addition of two capsules or the exact translation of a capsule by a 
%    vector
%
% Syntax:  
%    C = plus(summand1,summand2)
%
% Inputs:
%    summand1 - capsule or numerical vector
%    summand2 - capsule or numerical vector
%
% Outputs:
%    C - capsule after Minkowski addition
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    summand1 = C;
%    summand2 = [2; 2; 1];
%    C1 = C + summand1;
%    C2 = C + summand2;
%    figure; hold on;
%    plot(C); plot(C1); plot(C2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      04-March-2019 
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

%Find a capsule object
[C,summand] = findClassArg(summand1,summand2,'capsule');

try

    %Is summand a capsule?
    if isa(summand,'capsule')

        % add centers
        C.c = C.c + summand.c;
        % only replace generator if the one of the summand is longer
        length_C = norm(C.g);
        length_summand = norm(summand.g);
        if length_C < length_summand
            C.g = summand.g;
            radiusOfGenerator = length_C;
        else
            radiusOfGenerator = length_summand; 
        end
        % obtain new radius
        C.r = C.r + summand.r + radiusOfGenerator;
        
    %is summand a vector?
    elseif isnumeric(summand) && isvector(summand)
        %Calculate Minkowski sum
        C.c = C.c + summand;
        
    elseif isa(summand,'conPolyZono')
        
        C = summand + C;
    
    %something else?
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',summand1,summand2));
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if isempty(C)
        return
    elseif isemptyobject(summand)
        C = conZonotope(); return
    end

    % check whether different dimension of ambient space
    equalDimCheck(C,summand);

    % other error...
    rethrow(ME);

end

%------------- END OF CODE --------------