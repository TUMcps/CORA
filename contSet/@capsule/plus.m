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

% Authors:       Matthias Althoff
% Written:       04-March-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
        
    elseif isnumeric(summand) && isvector(summand)
        %Calculate Minkowski sum
        C.c = C.c + summand;
        
    elseif isa(summand,'conPolyZono')
        
        C = summand + C;

    elseif isa(summand,'interval')
        % shift center
        c_ = C.c + center(summand);
        % enlarge radius by radius of enclosing hyperball of interval
        r_ = C.r + radius(summand);
        C = capsule(c_,C.g,r_);

    elseif isa(summand,'zonotope')
        % outer-approximate zonotope by interval and use interval method
        C = C + interval(summand);
    
    elseif representsa_(summand,'origin',eps)
        % addition of only the origin (note: for some set representations,
        % this is a bit slow, so we put it at the end of the list...)
        return

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
    if representsa_(C,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',eps)
        C = capsule.empty(dim(C)); return
    end

    % check whether different dimension of ambient space
    equalDimCheck(C,summand);

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
