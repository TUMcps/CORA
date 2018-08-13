function [P] = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
% a mptPolytope
%
% Syntax:  
%    [P] = mtimes(matrix,P)
%
% Inputs:
%    factor1 - numerical matrix/interval matrix/mptPolytope object
%    factor2 - numerical matrix/interval matrix/mptPolytope object
%
% Outputs:
%    P - mptPolytope object
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
% Written:      01-February-2011
% Last update:  10-September-2015
%               15-June-2016
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%Find a zonotope object
%Is factor1 a mptPolytope?
if strcmp('mptPolytope',class(factor1))
    %initialize resulting polytope
    P=factor1;
    %initialize other summand
    matrix=factor2;
%Is factor2 a mptPolytope?    
elseif strcmp('mptPolytope',class(factor2))
    %initialize resulting zonotope
    P=factor2;
    %initialize other summand
    matrix=factor1;  
end

%get dimension
dim = dimension(P);

%define vector of zeros
z = zeros(dim,1);

%numeric matrix
if isnumeric(matrix)
    %convert scalar to matrix if necessary
    if length(matrix)==1
        matrix = matrix*eye(dim);
    end
    try
        try %MPT3
            P.P = matrix*P.P;
        catch %MPT2
            P.P=range(P.P,matrix,z);
        end
    catch
        disp('mapping via vertices');
        V = extreme(P.P);
        Vnew = matrix*V';
        P.P = polytope(Vnew');
   end
    
    
else
    %interval matrix 
    if strcmp('interval',class(matrix))
        %get minimum and maximum
        M_min=infimum(matrix);
        M_max=supremum(matrix);
    elseif strcmp('intervalMatrix',class(matrix))
        %get minimum and maximum
        M_min=infimum(matrix.int);
        M_max=supremum(matrix.int);
    end
    %get center of interval matrix
    T=0.5*(M_max+M_min);
    %get symmetric interval matrix
    Sval = 0.5*(M_max-M_min);
    S = interval(-Sval,Sval);

    %compute interval of polytope
    I = interval(P);

    %polytope of interval computations
    Iadd = S*I;
    Padd = mptPolytope(Iadd);

    %compute new polytope
    P=T*P + Padd; 
end


%------------- END OF CODE --------------