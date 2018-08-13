function [p] = randPointExtreme(qZ)
% randPoint - generates an extreme point within a quadZonotope when order=1
%
% Syntax:  
%    [p] = randPointExtreme(qZ)
%
% Inputs:
%    qZ - quadZonotope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      06-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get number of dependent and independent generators 
depGens = length(qZ.G(1,:));
indGens = length(qZ.Grest(1,:));

%random vectors for generator factors
beta = sign(1-2*rand(depGens,1));
gamma = sign(1-2*rand(indGens,1));

%start at center
p = qZ.c;

%first generator set
for i = 1:length(beta)
    p = p + beta(i)*qZ.G(:,i);
end

%second generator set
for i = 1:length(beta)
    p = p + beta(i)^2*qZ.Gsquare(:,i);
end

%third generator set
counter = 0;
for i = 1:(length(beta) - 1)
    for j = (i+1):length(beta)
        counter = counter + 1;
        p = p + beta(i)*beta(j)*qZ.Gquad(:,counter);
    end
end

%fourth generator set
for i = 1:length(gamma)
    p = p + gamma(i)*qZ.Grest(:,i);
end



%------------- END OF CODE --------------