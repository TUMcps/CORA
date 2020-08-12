function probZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%        probabilistic zonotope with another summand according to [1,(4)]
%
% Syntax:  
%    probZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - probabilistic zonotope object or other summand
%    summand2 - probabilistic zonotope object or other summand
%
% Outputs:
%    pZ - probabilistic zonotope 
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ1 = probZonotope(Z1,Z2);
% 
%    Z1 = [8 -1 0; 1 -2 -1];
%    Z2 = [0.4 -1.2; 0.2 -0.8];
%    probZ2 = probZonotope(Z1,Z2);
% 
%    probZ1 + probZ2
%
% References:
%    [1] M. Althoff et al. "Safety assessment for stochastic linear systems 
%        using enclosing hulls of probability density functions", ECC 2009
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      06-September-2007
% Last update:  24-August-2016
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

%Find a probabilistic zonotope object
%Is summand1 a zonotope?
if isa(summand1,'probZonotope')
    %initialize resulting zonotope
    probZ=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a zonotope?    
elseif isa(summand2,'probZonotope')
    %initialize resulting zonotope
    probZ=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand a probabilistic zonotope?
if isa(summand,'probZonotope')
    %Calculate minkowski sum
    probZ.Z=[probZ.Z(:,1)+summand.Z(:,1),probZ.Z(:,2:end),summand.Z(:,2:end)];
    %pZ.g=[pZ.g,summand.g];
    %pZ.sigma=[pZ.sigma,summand.sigma];
    probZ.cov=probZ.cov+summand.cov;

%Is summand a zonotope?
elseif isa(summand,'zonotope')
    %Calculate minkowski sum
    summandZ=summand.Z;
    probZ.Z=[probZ.Z(:,1)+summandZ(:,1),probZ.Z(:,2:end),summandZ(:,2:end)];
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate minkowski sum
    probZ.Z(:,1)=probZ.Z(:,1)+summand;
    
%something else?    
else
    % throw error for given arguments
    error(noops(summand1,summand2));
end

%------------- END OF CODE --------------