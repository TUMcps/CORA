function Z = plus(summand1,summand2)
% plus - overloaded '+' operator for the Minkowski addition of two
%    zonotopes or a zonotope with a vector
%
% Syntax:  
%    Z = plus(summand1,summand2)
%
% Inputs:
%    summand1 - zonotope object or numerical vector
%    summand2 - zonotope object or numerical vector
%
% Outputs:
%    Z - zonotope after Minkowski addition
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    summand1=Z;
%    summand2=[2; 2];
%    Z1=Z+summand1;
%    Z2=Z+summand2;
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'g');
%
% References:
%    [1] M. Althoff. "Reachability analysis and its application to the 
%        safety assessment of autonomous cars", 2010
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
%               04-March-2019
%               13-August-2019
% Last revision:---

%------------- BEGIN CODE --------------

    % determine zonotope object
    if isa(summand1,'zonotope')
        Z = summand1;
        summand = summand2;  
    elseif isa(summand2,'zonotope')
        Z = summand2;
        summand = summand1;  
    end

    % different cases depending on the class of the second summand
    if isa(summand,'zonotope')

        % see Equation 2.1 in [1]
        Z.Z(:,1) = Z.Z(:,1)+summand.Z(:,1);
        Z.Z(:,(end+1):(end+length(summand.Z(1,2:end)))) = summand.Z(:,2:end);

    elseif isnumeric(summand)

        Z.Z(:,1) = Z.Z(:,1)+summand;

    elseif isa(summand,'interval')

        Z = Z + zonotope(summand);

    elseif isa(summand,'mptPolytope') || isa(summand,'conZonotope') || ...
           isa(summand,'zonoBundle') || isa(summand,'polyZonotope') || ...
           isa(summand,'conPolyZono')

        Z = summand + Z;        

    else
        % throw error for given arguments
        error(noops(summand1,summand2));
    end

%------------- END OF CODE --------------