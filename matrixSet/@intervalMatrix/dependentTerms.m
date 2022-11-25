function [intSq,intH] = dependentTerms(obj,r)
% dependentTerms - computes exact Taylor terms of an interval matrix square
%    and an interval matrix exponential
%
% These different tasks are computed in one m-file to save computation
% time: the for loop has to be executed only once and help functions do not
% have to be called so often
%
% Syntax:  
%    [intSq,intH] = dependentTerms(obj,r)
%
% Inputs:
%    obj - linear interval system (linIntSys) object
%    r - time step increment
%
% Outputs:
%    intSq - exact square matrix
%    intH - exact Taylor terms up to second order
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Matthias Althoff
% Written:      20-February-2007 
% Last update:  30-April-2007
%               23-September-2010
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object structure
A=obj.int;
dim=obj.dim;

%initialize the square of A (sq), the first term of the interval
%exponential (H)
sq=0*A;
H=0*A;
I=eye(dim); %identity matrix

%get diagonal elements of A (diagA)
diagA=interval(zeros(dim),zeros(dim)); %initialize
for i=1:dim
    diagA(i,i)=A(i,i);
end

%compute elements of sq and H
for i=1:dim
    %i neq j
    %auxiliary value s
    s=sum(A,i);
    %auxiliary value b
    b=A(i,:); b(i)=0;
    %auxiliary matrix C
    C=I*A(i,i)+diagA;
    %compute non-diagonal elements of sq
    sq(i,:)=b*C+s;
    %compute non-diagonal elements of H
    H(i,:)=b*(I*r+0.5*C*r^2)+0.5*r^2*s;
    
    %i=j
    %compute diagonal elements of sq
    sq(i,i)=sq(i,i)+A(i,i)^2;            
    %auxiliary values for H, Hu
    a_inf=infimum(A(i,i));
    a_sup=supremum(A(i,i));
    %compute diagonal elements for H
    kappa=max(a_inf*r+0.5*a_inf^2*r^2, a_sup*r+0.5*a_sup^2*r^2);
    H(i,i)=H(i,i)+interval(g(A(i,i),r),kappa);         
    %--------------------------------------------------------------
end

%write as interval matrices
intSq=intervalMatrix(r^2*center(sq),r^2*rad(sq));
intH=intervalMatrix(center(H)+I,rad(H));

end


% AUXILIARY FUNCTIONS

function [res]=g(a,r)
    if isIntersecting(interval(-1/r,-1/r),a)
        res=-0.5;
    else
        a_inf=infimum(a);
        a_sup=supremum(a);
        res=min(a_inf*r+0.5*a_inf^2*r^2, a_sup*r+0.5*a_sup^2*r^2);
    end
end


function [res]=gu(a,r)
    if in(interval(-3/(2*r)),a)
        res=-9/24*r;
    else
        a_inf=infimum(a);
        a_sup=supremum(a);
        res=min(0.5*a_inf*r^2+1/6*a_inf^2*r^3,...
            0.5*a_sup*r^2+1/6*a_sup^2*r^3);
    end    
end

    
function s=sum(A,i)
%sum function:
%s=0.5 \sum_{k:k\neq i,k\neq j} a_{ik}a_{kj}t^2
    n=length(A);
    k=0:n;
    ind=k*n+1:n+1:n^2;
    A(ind)=zeros(n,1);

    s=A(i,:)*A;
end

%------------- END OF CODE --------------