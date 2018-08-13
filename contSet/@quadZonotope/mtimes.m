function [Z] = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
% interval matrix with a quadZonotope
%
% Syntax:  
%    [Z] = mtimes(matrix,Z)
%
% Inputs:
%    matrix - numerical or interval matrix
%    Z - quadZonotope object 
%
% Outputs:
%    Z - quadZonotpe after multiplication of a matrix with a quadZonotope
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      04-September-2012 
% Last update:  18-March-2016
%               19-November-2017
% Last revision:---

%------------- BEGIN CODE --------------

%Find a zonotope object
%Is factor1 a zonotope?
if isa(factor1,'quadZonotope')
    %initialize resulting zonotope
    Z=factor1;
    %initialize other summand
    matrix=factor2;
%Is factor2 a zonotope?    
elseif isa(factor2,'quadZonotope')
    %initialize resulting zonotope
    Z=factor2;
    %initialize other summand
    matrix=factor1;  
end

%numeric matrix
if isnumeric(matrix)
    Z.c=matrix*Z.c;
    if ~isempty(Z.G)
        Z.G=matrix*Z.G;
    end
    if ~isempty(Z.Gsquare)
        Z.Gsquare=matrix*Z.Gsquare;
    end
    if ~isempty(Z.Gquad)
        Z.Gquad=matrix*Z.Gquad;
    end
    if ~isempty(Z.Grest)
        Z.Grest=matrix*Z.Grest;
    end

 
%interval matrix
elseif isa(matrix,'interval')
    %get minimum and maximum
    M_min=infimum(matrix);
    M_max=supremum(matrix);
    %get center of interval matrix
    T=0.5*(M_max+M_min);
    %get symmetric interval matrix
    S=0.5*(M_max-M_min);
    
    %absolute sum of all generators
    Zabssum=sum(abs([Z.c, Z.G, Z.Gsquare, Z.Gquad, Z.Grest]), 2);
    
    %compute new zonotope
    Z.c = T*Z.c;
    Z.G = T*Z.G;
    Z.Gsquare = T*Z.Gsquare;
    Z.Gquad = T*Z.Gquad;
    Z.Grest = [T*Z.Grest, diag(S*Zabssum)];
    
    
%interval matrix 
elseif isa(matrix,'intervalMatrix')
    %get minimum and maximum
    M_min=infimum(matrix.int);
    M_max=supremum(matrix.int); 
    %get center of interval matrix
    T=0.5*(M_max+M_min);
    %get symmetric interval matrix
    S=0.5*(M_max-M_min);
    
    %absolute sum of all generators
    Zabssum=sum(abs([Z.c, Z.G, Z.Gsquare, Z.Gquad, Z.Grest]), 2);
    
    %compute new zonotope
    Z.c = T*Z.c;
    Z.G = T*Z.G;
    Z.Gsquare = T*Z.Gsquare;
    Z.Gquad = T*Z.Gquad;
    Z.Grest = [T*Z.Grest, diag(S*Zabssum)];

    
%matrix zonotope 
elseif isa(matrix,'matZonotope')
    
    disp('not yet implemented')
end    




%------------- END OF CODE --------------