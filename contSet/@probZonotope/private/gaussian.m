function f=gaussian(x,Sigma)
% gaussian - computes the values of a Gaussian (normal) distribution
%
% Syntax:
%    f=gaussian(x,Sigma)
%
% Inputs:
%    x -
%    Sigma -
%
% Outputs:
%    f - 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       10-October-2007
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%get dimension
d=length(Sigma);

%acceleration for 2-dim plots
if d==2
    x1=x(1,:);
    x2=x(2,:);
    auxVar=1/det(Sigma)*(Sigma(2,2)*x1.^2-2*Sigma(1,2)*x1.*x2+Sigma(1,1)*x2.^2);
    f=1/((2*pi)^(d/2)*det(Sigma)^(1/2))*exp(-1/2*auxVar);
%otherwise
else
    f=1/((2*pi)^(d/2)*det(Sigma)^(1/2))*exp(-1/2*x'*inv(Sigma)*x);
end

% ------------------------------ END OF CODE ------------------------------
