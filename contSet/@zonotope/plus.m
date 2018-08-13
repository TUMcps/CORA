function [Z] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of two
% zonotopes or a zonotope with a vector
%
% Syntax:  
%    [Z] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - zonotope object or numerical vector
%    summand2 - zonotope object or numerical vector
%
% Outputs:
%    Z - Zonotpe after Minkowsi addition
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    summand1=Z;
%    summand2=[2; 2];
%    Z1=Z+summand1;
%    Z2=Z+summand2;
%    plot(Z);
%    hold on
%    plot(Z1);
%    plot(Z2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      30-September-2006 
% Last update:  23-March-2007
%               14-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%Find a zonotope object
%Is summand1 a zonotope?
if strcmp('zonotope',class(summand1))
    %initialize resulting zonotope
    Z=summand1;
    %initialize other summand
    summand=summand2;
%Is summand2 a zonotope?    
elseif strcmp('zonotope',class(summand2))
    %initialize resulting zonotope
    Z=summand2;
    %initialize other summand
    summand=summand1;  
end

%Is summand a zonotope?
if strcmp('zonotope',class(summand))
    %Calculate minkowski sum
    Z.Z(:,1)=Z.Z(:,1)+summand.Z(:,1);
    Z.Z(:,(end+1):(end+length(summand.Z(1,2:end)))) = summand.Z(:,2:end);
    %Z.Z=[Z.Z(:,1)+summand.Z(:,1),Z.Z(:,2:end),summand.Z(:,2:end)];
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate minkowski sum
    Z.Z(:,1)=Z.Z(:,1)+summand;
    
%something else?    
else
    Z.Z=[];
    disp('this operation is not implemented');
end

%------------- END OF CODE --------------