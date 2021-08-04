function cZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%        constrained zonotope with other set representations
%
% Syntax:  
%    cZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - conZonotope object or numerical vector
%    summand2 - conZonotope object or numerical vector
%
% Outputs:
%    cZ - constrained Zonotope after Minkowski addition
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZono = conZonotope(Z,A,b);
%    cPlus = cZono + [5;4];
%
%    hold on
%    plot(cZono,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(cPlus,[1,2],'b','Filled',true,'EdgeColor','none');
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
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% find a conZonotope object
if isa(summand1, 'conZonotope')
    cZ=summand1;
    summand=summand2;
elseif isa(summand2, 'conZonotope')
    cZ=summand2;
    summand=summand1;  
end

% handle different classes of the second summand
if isa(summand, 'conZonotope')
    
    % Calculate minkowski sum (Equation (12) in reference paper [1])
    cZ.Z(:,1)=cZ.Z(:,1)+summand.Z(:,1);
    cZ.Z(:,(end+1):(end+length(summand.Z(1,2:end)))) = summand.Z(:,2:end);
    
    if isempty(cZ.A)
        if ~isempty(summand.A)
            cZ.A = [zeros(size(summand.A,1), ...
                   size(cZ.Z,2)-1-size(summand.A,2)),summand.A]; 
            cZ.b = summand.b;
        end
    else
        if isempty(summand.A)
            cZ.A = [cZ.A,zeros(size(cZ.A,1),size(summand.Z,2)-1)];
        else
            cZ.A = blkdiag(cZ.A, summand.A);
            cZ.b = [cZ.b; summand.b];
        end
    end
    
    cZ.ksi = [];
    cZ.R = [];
    
elseif isnumeric(summand) 
    
    cZ.Z(:,1) = cZ.Z(:,1) + summand;
    
elseif isa(summand,'zonotope') || isa(summand,'interval') || ...
       isa(summand,'mptPolytope') || isa(summand,'zonoBundle')
    
    cZ = cZ + conZonotope(summand);
    
elseif isa(summand,'polyZonotope') || isa(summand,'conPolyZono')
    
    cZ = polyZonotope(cZ) + summand;
     
else
    % throw error for given arguments
    error(noops(summand1,summand2));
end

%------------- END OF CODE --------------