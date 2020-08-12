function fSum=mixedGaussian(x,omega,Sigma,mu)
% mixedGaussian - 
%
% Syntax:  
%    fSum=mixedGaussian(x,omega,Sigma,mu)
%
% Inputs:
%    x -
%    omega - 
%    Sigma -
%    mu - 
%
% Outputs:
%    fSum - 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      ---
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

d=length(x);

for i=1:length(Sigma)
    f(i,:)=omega{i}/((2*pi)^(d/2)*det(Sigma{i})^(1/2))*exp(-1/2*(x-mu{i})'*inv(Sigma{i})*(x-mu{i}));
end
fSum=sum(f);

%------------- END OF CODE --------------