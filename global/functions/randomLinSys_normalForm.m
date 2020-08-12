function [linSys]=randomLinSys_normalForm(dim,realInterval,imagInterval)
% randomLinSys - generates a random linear system by randomly placing poles
%
% Syntax:  
%    [linSys]=randomLinSys(dim,realInterval,imagInterval,type)
%
% Inputs:
%    dim - dimension
%    realInterval - interval of real values of poles
%    imagInterval - interval of imag values of poles 
%    type - type of probability distribution
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


%create system matrix A
A = zeros(dim);
for i=1:(dim-1)
    for j=1:i
        A(dim-i,dim+1-j) = 0.03*(1-2*rand(1,1));
    end
end
for i=1:(dim)
    A(i,i) = -3*rand(1,1)-0.2;
end

B = 2*rand(dim)-1;

%instantiate linear system
linSys = linearSys('linSys',A,B);

%------------- END OF CODE --------------