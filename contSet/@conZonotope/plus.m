function [Z] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%        constrained zonotope with other set representations
%
% Syntax:  
%    [Z] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - conZonotope object or numerical vector
%    summand2 - conZonotope object or numerical vector
%
% Outputs:
%    Z - constrained Zonotope after Minkowsi addition
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZono = conZonotope(Z,A,b);
%    cPlus = cZono + [5;4];
%
%    hold on
%    plotFilled(cZono,[1,2],'r','EdgeColor','none');
%    plotFilled(cPlus,[1,2],'b','EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      05-December-2017 
% Last update:  15-May-2018
% Last revision:---

%------------- BEGIN CODE --------------

% Find a conZonotope object
if isa(summand1, 'conZonotope')
    Z=summand1;
    summand=summand2;
elseif isa(summand2, 'conZonotope')
    Z=summand2;
    summand=summand1;  
end

% Handle different classes of the second summand
if isa(summand, 'conZonotope')
    
    % Calculate minkowski sum (Equation (12) in reference paper [1])
    Z.Z(:,1)=Z.Z(:,1)+summand.Z(:,1);
    Z.Z(:,(end+1):(end+length(summand.Z(1,2:end)))) = summand.Z(:,2:end);
    Z.A = blkdiag(Z.A, summand.A);
    Z.b = [Z.b; summand.b];
    
    Z.ksi = [];
    Z.R = [];
    
elseif isa(summand, 'zonotope')
    
    % Calculate minkowski sum (Equation (12) in reference paper [1])
    Z.Z(:,1)=Z.Z(:,1)+summand.Z(:,1);
    Z.Z(:,(end+1):(end+length(summand.Z(1,2:end)))) = summand.Z(:,2:end);
    Z.A = [Z.A, zeros(size(Z.A,1),size(summand.Z,2)-1)];

    Z.ksi = [];
    Z.R = [];
    
elseif isa(summand, 'interval')
    
    % Convert interval to zontope 
    Z = Z + zonotope(summand);
    
elseif isnumeric(summand)       % summand is a vector
    
    %Calculate minkowski sum
    Z.Z(:,1)=Z.Z(:,1)+summand;
     
else
    error('This operation is not implemented');
end

%------------- END OF CODE --------------