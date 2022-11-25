function Asq = exactSquare(Aint)
% exactSquare - computes the exact square of an interval matrix
%
% The computation is according to [1].
% 
% A more easily understandable implementation can be found in the unit
% test 'test_intervalMatrix_exactSquare'.
%
% Syntax:  
%    Asq = exactSquare(Aint)
%
% Inputs:
%    A - interval matrix 
%
% Outputs:
%    Asq - resulting interval matrix
%
% Example:
%    -
%
% References:
%   [1] Kosheleva, O.; Kreinovich, V.; Mayer, G. & Nguyen, H. T. Computing the 
% cube of an interval matrix is NP-Hard Proc. of the ACM symposium on 
% Applied computing, 2005, 1449-1453
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 
%
% Author:       Matthias Althoff
% Written:      04-January-2009 
% Last update:  02-November-2017
% Last revision:---

%------------- BEGIN CODE --------------

%extract intervals
A = Aint.int;

%obtain dimension
dim=length(A);

%initialize
sq=0*A;
E=eye(dim); %identity matrix

%compute result for diagonal and non-diagonal elements
%compute elements of H, Hu and sq
for i=1:dim
    %i neq j
    %auxiliary value s
    s=sum(A,i);
    %auxiliary value b
    b=A(i,:); b(i)=0;
    %auxiliary matrix C
    C=E*A(i,i)+diag(diag(A));
    %compute non-diagonal elements of sq
    sq(i,:)=b*C+s;
    
    %i=j
    %compute diagonal elements of sq
    sq(i,i)=sq(i,i)+A(i,i)^2;            
end

% include result in intervalMatrix object
Asq = intervalMatrix([]);
Asq.int = sq;

%sum function:
%s=\sum_{k:k\neq i,k\neq j} a_{ik}a_{kj}
function s=sum(A,i)

% for k=1:length(A)
%     A(k,k)=0;
% end

%get indices that should be 0
n=length(A);
k=0:n;
ind=k*n+1:n+1:n^2;
A(ind)=zeros(n,1);

s=A(i,:)*A;

%------------- END OF CODE --------------