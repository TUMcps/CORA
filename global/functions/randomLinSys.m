function [linSys]=randomLinSys(dim,realInterval,imagInterval)
% randomLinSys - generates a random linear system by randomly placing poles
%
% Syntax:  
%    [linSys] = randomLinSys(dim,realInterval,imagInterval)
%
% Inputs:
%    dim - dimension
%    realInterval - array(2) lower + upper bound of real values of poles
%    imagInterval - array(2) lower + upper bound of imag values of poles 
%
% Outputs:
%    linSys - linear system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      11-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%determine number of real and conjugate complex eigenvalues
if rem(dim,2)>0
    nReal = 1;
else
    nReal = 0;
end
nConj = floor(dim/2);


%create real values
for i=1:(nReal+nConj)
    %uniform distribution
    realVal(i) = unifrnd(realInterval(1),realInterval(2));
end
%create imag values
for i=1:nConj
    %uniform distribution
    imagVal(i) = unifrnd(imagInterval(1),imagInterval(2));
end


%generate pole vector
for i=1:nConj
    p(2*i-1) = complex(realVal(i),imagVal(i));
    p(2*i) = complex(realVal(i),-imagVal(i));
end
for i=1:nReal
    %uniform distribution
    p(end+1) = realVal(end);
end

%create system matrix A
A = zeros(dim);
for i=1:(dim-1)
    A(i,i+1) = 1;
end
coeff = poly(p);
ind = linspace(dim+1,2,dim);
A(dim,:) = -coeff(ind);

%obtain orthogonal random matrix to transform system
V = 2*rand(dim)-1;
%orthogonalize
[Q,R] = qr(V);

%Q = gallery('randhess',dim);

%create random A, B
A = Q'*A*Q;
B = 2*rand(dim)-1;

%instantiate linear system
linSys = linearSys('linSys',A,B);

%------------- END OF CODE --------------