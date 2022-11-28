function Z = mtimes(factor1,factor2)
% mtimes - overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a zonotope
%
% Syntax:  
%    Z = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical or interval matrix
%    factor2 - zonotope object 
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z = zonotope([1 1 0; 0 0 1]);
%    M = [0 1; 1 0];
%    Zmat = M*Z;
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(Z,[1,2],'r');
%
% References:
%    [1] M. Althoff. "Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010
%    [2] M. Althoff et al. "Modeling, Design, and Simulation of Systems 
%        with Uncertainties". 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      30-September-2006 
% Last update:  07-September-2007
%               05-January-2009
%               06-August-2010
%               01-February-2011
%               08-February-2011
%               18-November-2015
% Last revision:---

%------------- BEGIN CODE --------------

%Find a zonotope object
[Z,M] = findClassArg(factor1,factor2,'zonotope');

try

    %numeric matrix
    if isnumeric(M)
        Z.Z=M*Z.Z;
        
    % interval (see Theorem 3.3 in [1])
    elseif isa(M,'interval')
        %get minimum and maximum
        M_min=infimum(M);
        M_max=supremum(M);
        %get center of interval matrix
        T=0.5*(M_max+M_min);
        %get symmetric interval matrix
        S=0.5*(M_max-M_min);
        Zabssum=sum(abs(Z.Z),2);
        %compute new zonotope
        Z.Z=[T*Z.Z,diag(S*Zabssum)]; 
    
    % interval matrix (see Theorem 3.3 in [1])
    elseif isa(M,'intervalMatrix')
        %get minimum and maximum
        M_min=infimum(M.int);
        M_max=supremum(M.int); 
        %get center of interval matrix
        T=0.5*(M_max+M_min);
        %get symmetric interval matrix
        S=0.5*(M_max-M_min);
        Zabssum=sum(abs(Z.Z),2);
        %compute new zonotope
        Z.Z=[T*Z.Z,diag(S*Zabssum)]; 
    
    % matrix zonotope (see Sec. 4.4.1 in [2])
    elseif isa(M,'matZonotope')
        %obtain first zonotope
        Znew=M.center*Z.Z;
        %compute further zonotopes and add them up
        for i=1:M.gens
            Zadd=M.generator{i}*Z.Z;
            Znew(:,(end+1):(end+length(Z.Z(1,:))))=Zadd;
        end
        %write to Z.Z
        Z.Z=Znew;
    
    else
        % throw error
        throw(CORAerror('CORA:noops',M,Z));
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
    end

    % check whether different dimension of ambient space
    equalDimCheck(Z,M);

    % other error...
    rethrow(ME);

end

%------------- END OF CODE --------------