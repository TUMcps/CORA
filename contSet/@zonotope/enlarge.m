function Z = enlarge(Z,factorVec)
% enlarge - enlarges the generators of a zonotope by a vector of factors
%
% Syntax:  
%    Z = enlarge(Z,factorVec)
%
% Inputs:
%    Z - zonotope object
%    factorVec - column vector of factors for the enlargement of each dimension
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z = zonotope([0;0],rand(2,5));
%    Zlarger = enlarge(Z,1.2);
%    Zsmaller = enlarge(Z,[0.9;0.8]);
% 
%    plot(Z); hold on;
%    plot(Zlarger,[1,2],'g');
%    plot(Zsmaller,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      20-November-2010 
% Last update:  21-April-2020 (.* instead of diag())
% Last revision:---

%------------- BEGIN CODE --------------

Z.Z(:,2:end) = factorVec .* Z.Z(:,2:end);

%------------- END OF CODE --------------