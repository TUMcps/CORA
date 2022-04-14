function cPZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%        constrained polynomial zonotope with other set representations
%
% Syntax:  
%    cPZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - conPolyZono object, other set rep., or numerical vector
%    summand2 - conPolyZono object, other set rep., or numerical vector
%
% Outputs:
%    cPZ - conPolyZonotope after Minkowsi addition
%
% Example: 
%    c = [0;0];
%    G = [2 2;2 -1];
%    expMat = [1 0; 0 1];
%    A = [1 1];
%    b = 0;
%    expMat_ = [2 0; 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    E = ellipsoid(0.1*[2 1;1 2]);
%
%    res = cPZ + E; 
%
%    figure; hold on;
%    plot(res,[1,2],'b','Filled',true,'EdgeColor','none','Splits',25);
%    plot(E,[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, exactPlus

% Author:       Niklas Kochdumper
% Written:      14-August-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % determine which summand is the conPolyZono object
    if isa(summand1,'conPolyZono')
        cPZ = summand1;
        summand = summand2;  
    elseif isa(summand2,'conPolyZono')
        cPZ = summand2;
        summand = summand1;  
    end
    
    % different cases for different set represnetations
    if isnumeric(summand) || isa(summand,'interval') || ...
       isa(summand,'zonotope')
    
       Z = zonotope(summand);
       cPZ.c = cPZ.c + center(Z);
       cPZ.Grest = [cPZ.Grest,generators(Z)];
       
    elseif isa(summand,'conPolyZono') || isa(summand,'mptPolytope') || ...
           isa(summand,'zonoBundle') || isa(summand,'conZonotope') || ...
           isa(summand,'ellipsoid') || isa(summand,'capsule') || ...
           isa(summand,'polyZonotope') || isa(summand,'taylm')
        
         % convert to conPolyZono object
         summand = conPolyZono(summand);
       
         % update states
         cPZ.c = cPZ.c + summand.c;
         cPZ.G = [cPZ.G, summand.G];
         cPZ.expMat = blkdiag(cPZ.expMat,summand.expMat);
         
         cPZ.Grest = [cPZ.Grest, summand.Grest];
         
         % update constraints
         cPZ = updateConstraints(cPZ,cPZ,summand);
         
    else
         % throw error for given arguments
         error(noops(cPZ,summand));
    end
end

%------------- END OF CODE --------------