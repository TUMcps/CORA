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

% Authors:       Matthias Althoff
% Written:       30-September-2006 
% Last update:   23-March-2007
%                14-August-2016
%                04-March-2019
%                13-August-2019
%                14-February-2024 (MW, prevent sum with row vectors/matrices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine zonotope object
[Z,summand] = findClassArg(summand1,summand2,'zonotope');

try

    % different cases depending on the class of the second summand
    if isa(summand,'zonotope')
    
        % see Equation 2.1 in [1]
        Z.c = Z.c+summand.c;
        if isempty(Z.c)
            Z.G = zeros(dim(Z),0);
        else
            Z.G(:,(end+1):(end+size(summand.G,2))) = summand.G;
        end
    
    elseif isnumeric(summand) && (isscalar(summand) || all(size(summand) == [dim(Z),1]))
        % summand has to be a scalar or a column vector of correct size
        
        Z.c = Z.c + summand;
    
    elseif isa(summand,'interval')
    
        Z = Z + zonotope(summand);

    elseif isa(summand,'polytope') || isa(summand,'conZonotope') || ...
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

    % check whether different dimension of ambient space
    equalDimCheck(Z,summand);

    % check for empty sets
    if representsa_(Z,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',1e-10)
        Z = zonotope.empty(dim(Z)); return
    end

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
