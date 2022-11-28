function Z = plus(summand1,summand2)
% plus - overloaded '+' operator for the Minkowski addition of two
%    zonotopes or a zonotope with a vector
%
% Syntax:  
%    Z = plus(summand1,summand2)
%
% Inputs:
%    summand1 - zonotope object or numerical vector
%    summand2 - zonotope object, contSet object, or numerical vector
%
% Outputs:
%    Z - zonotope after Minkowski addition
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    summand1=Z;
%    summand2=[2; 2];
%    Z1=Z+summand1;
%    Z2=Z+summand2;
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'g');
%
% References:
%    [1] M. Althoff. "Reachability analysis and its application to the 
%        safety assessment of autonomous cars", 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      30-September-2006 
% Last update:  23-March-2007
%               14-August-2016
%               04-March-2019
%               13-August-2019
% Last revision:---

%------------- BEGIN CODE --------------

% determine zonotope object
[Z,summand] = findClassArg(summand1,summand2,'zonotope');

try

    % different cases depending on the class of the second summand
    if isa(summand,'zonotope')
    
        % see Equation 2.1 in [1]
        Z.Z(:,1) = Z.Z(:,1)+summand.Z(:,1);
        Z.Z(:,(end+1):(end+length(summand.Z(1,2:end)))) = summand.Z(:,2:end);
    
    elseif isnumeric(summand)
    
        Z.Z(:,1) = Z.Z(:,1)+summand;
    
    elseif isa(summand,'interval')
    
        Z = Z + zonotope(summand);
    
    elseif isa(summand,'mptPolytope') || isa(summand,'conZonotope') || ...
           isa(summand,'zonoBundle') || isa(summand,'polyZonotope') || ...
           isa(summand,'conPolyZono')
    
        Z = summand + Z;        
    
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
    if isempty(Z)
        return
    elseif (isnumeric(summand) && isempty(summand)) || ...
            (isa(summand,'contSet') && isemptyobject(summand))
        Z = zonotope(); return
    end

    % check whether different dimension of ambient space
    equalDimCheck(Z,summand);

    % other error...
    rethrow(ME);

end

%------------- END OF CODE --------------