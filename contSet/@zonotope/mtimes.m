function [Z] = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
% interval matrix with a zonotope
%
% Syntax:  
%    [Z] = mtimes(matrix,Z)
%
% Inputs:
%    matrix - numerical or interval matrix
%    Z - zonotope object 
%
% Outputs:
%    Z - Zonotpe after multiplication of a matrix with a zonotope
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    matrix=[0 1; 1 0];
%    plot(Z);
%    hold on
%    Z=matrix*Z;
%    plot(Z);
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
    %update orientation
    Z = updateOrientation(Z,matrix);
    
    
%interval matrix
elseif strcmp('interval',class(matrix))
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
    %update orientation
    Z = updateOrientation(Z,T);

    
%interval matrix 
elseif strcmp('intervalMatrix',class(matrix))
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
    %update orientation
    Z = updateOrientation(Z,T);

    
%matrix zonotope 
elseif strcmp('matZonotope',class(matrix))
    %obtain first zonotope
    Znew=matrix.center*Z.Z;
    %compute further zonotopes and add them up
    for i=1:matrix.gens
        Zadd=matrix.generator{i}*Z.Z;
        Znew(:,(end+1):(end+length(Z.Z(1,:))))=Zadd;
    end
    %write to Z.Z
    Z.Z=Znew;
%     %update orientation
%     Z = updateOrientation(Z,matrix.center);
end    
% %matrix zonotope 
% elseif strcmp('matZonotope',class(matrix))
%     %preallocate
%     rows = length(matrix.center(:,1));
%     cols_add = length(Z.Z(1,:));
%     cols = (matrix.gens+1)*cols_add;
%     Znew = zeros(rows,cols);
%     %obtain first zonotope
%     Znew(:,1:cols_add)=matrix.center*Z.Z;
%     %compute further zonotopes and add them up
%     for i=1:matrix.gens
%         Znew(:,(i*cols_add+1):((i+1)*cols_add))=matrix.generator{i}*Z.Z;
%     end
%     %write to Z.Z
%     Z.Z=Znew;
%     %update orientation
%     %Z = updateOrientation(Z,matrix.center);
% end


function Z = updateOrientation(Z,M)

if ~isempty(Z.O)
    %linear map
    Z.O = M*Z.O;

    %normalization
    for i=1:length(M)
        Z.O(:,i) = Z.O(:,i)/norm(Z.O(:,i));
    end
end



%------------- END OF CODE --------------