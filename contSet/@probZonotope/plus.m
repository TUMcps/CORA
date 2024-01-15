function probZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%        probabilistic zonotope with another summand according to [1,(4)]
%
% Syntax:
%    probZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - probabilistic zonotope object or other summand
%    summand2 - probabilistic zonotope object or other summand
%
% Outputs:
%    pZ - probabilistic zonotope 
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ1 = probZonotope(Z1,Z2);
% 
%    Z1 = [8 -1 0; 1 -2 -1];
%    Z2 = [0.4 -1.2; 0.2 -0.8];
%    probZ2 = probZonotope(Z1,Z2);
% 
%    probZ1 + probZ2
%
% References:
%    [1] M. Althoff et al. "Safety assessment for stochastic linear systems 
%        using enclosing hulls of probability density functions", ECC 2009
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       06-September-2007
% Last update:   24-August-2016
%                05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine probabilistic zonotope object
[probZ,summand] = findClassArg(summand1,summand2,'probZonotope');
    
try

    %Is summand a probabilistic zonotope?
    if isa(summand,'probZonotope')
        %Calculate minkowski sum
        probZ.Z=[probZ.Z(:,1)+summand.Z(:,1),probZ.Z(:,2:end),summand.Z(:,2:end)];
        %pZ.g=[pZ.g,summand.g];
        %pZ.sigma=[pZ.sigma,summand.sigma];
        probZ.cov=probZ.cov+summand.cov;
    
    %Is summand a zonotope?
    elseif isa(summand,'zonotope')
        %Calculate minkowski sum
        summandZ=[summand.c,summand.G];
        probZ.Z=[probZ.Z(:,1)+summandZ(:,1),probZ.Z(:,2:end),summandZ(:,2:end)];
        
    %is summand a vector?
    elseif isnumeric(summand)
        %Calculate minkowski sum
        probZ.Z(:,1)=probZ.Z(:,1)+summand;
        
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

    % check whether different dimension of ambient space
    equalDimCheck(probZ,summand);

    % check for empty sets
    if representsa_(probZ,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',eps)
        probZ = probZonotope(); return
    end

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
