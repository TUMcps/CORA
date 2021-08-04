function Z = mtimes(factor1,factor2)
% mtimes - overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a zonotope
%
% Syntax:  
%    Z = mtimes(matrix,Z)
%
% Inputs:
%    matrix - numerical or interval matrix
%    Z - zonotope object 
%
% Outputs:
%    Z - zonotope after multiplication of a matrix with a zonotope
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    matrix=[0 1; 1 0];
%    Zmat = matrix*Z;
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
%Is factor1 a zonotope?
if isa(factor1,'zonotope')
    %initialize resulting zonotope
    Z=factor1;
    %initialize other summand
    matrix=factor2;
%Is factor2 a zonotope?    
elseif isa(factor2,'zonotope')
    %initialize resulting zonotope
    Z=factor2;
    %initialize other summand
    matrix=factor1;  
end

%numeric matrix
if isnumeric(matrix)
    Z.Z=matrix*Z.Z;
    
    
% interval (see Theorem 3.3 in [1])
elseif isa(matrix,'interval')
    %get minimum and maximum
    M_min=infimum(matrix);
    M_max=supremum(matrix);
    %get center of interval matrix
    T=0.5*(M_max+M_min);
    %get symmetric interval matrix
    S=0.5*(M_max-M_min);
    Zabssum=sum(abs(Z.Z),2);
    %compute new zonotope
    Z.Z=[T*Z.Z,diag(S*Zabssum)]; 

    
% interval matrix (see Theorem 3.3 in [1])
elseif isa(matrix,'intervalMatrix')
    %get minimum and maximum
    M_min=infimum(matrix.int);
    M_max=supremum(matrix.int); 
    %get center of interval matrix
    T=0.5*(M_max+M_min);
    %get symmetric interval matrix
    S=0.5*(M_max-M_min);
    Zabssum=sum(abs(Z.Z),2);
    %compute new zonotope
    Z.Z=[T*Z.Z,diag(S*Zabssum)]; 

    
% matrix zonotope (see Sec. 4.4.1 in [2])
elseif isa(matrix,'matZonotope')
    %obtain first zonotope
    Znew=matrix.center*Z.Z;
    %compute further zonotopes and add them up
    for i=1:matrix.gens
        Zadd=matrix.generator{i}*Z.Z;
        Znew(:,(end+1):(end+length(Z.Z(1,:))))=Zadd;
    end
    %write to Z.Z
    Z.Z=Znew;
end    



%------------- END OF CODE --------------