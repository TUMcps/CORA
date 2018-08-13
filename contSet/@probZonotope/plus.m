function [pZ] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
% probabilistic zonotope with another summand
%
% Syntax:  
%    [pZ] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - probabilistic zonotope object or other summand
%    summand2 - probabilistic zonotope object or other summand
%
% Outputs:
%    pZ - probabilistic Zonotpe 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      06-September-2007
% Last update:  24-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%Find a probabilistic zonotope object
%Is summand1 a zonotope?
if strcmp('probZonotope',class(summand1))
    %initialize resulting zonotope
    pZ=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a zonotope?    
elseif strcmp('probZonotope',class(summand2))
    %initialize resulting zonotope
    pZ=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand a probabilistic zonotope?
if strcmp('probZonotope',class(summand))
    %Calculate minkowski sum
    pZ.Z=[pZ.Z(:,1)+summand.Z(:,1),pZ.Z(:,2:end),summand.Z(:,2:end)];
    %pZ.g=[pZ.g,summand.g];
    %pZ.sigma=[pZ.sigma,summand.sigma];
    pZ.cov=pZ.cov+summand.cov;

%Is summand a zonotope?
elseif strcmp('zonotope',class(summand))
    %Calculate minkowski sum
    summandZ=get(summand,'Z');
    pZ.Z=[pZ.Z(:,1)+summandZ(:,1),pZ.Z(:,2:end),summandZ(:,2:end)];
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate minkowski sum
    pZ.Z(:,1)=pZ.Z(:,1)+summand;
    
%something else?    
else
    pZ.Z=[];
    disp('this operation is not implemented');
end

%------------- END OF CODE --------------