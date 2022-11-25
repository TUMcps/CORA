function cPZ = mtimes(matrix,cPZ)
% mtimes - Overloaded '.*' operator for the multiplication of a matrix with 
%          a constrained polynomial zonotope
%
% Syntax:  
%    cPZ = mtimes(matrix,cPZ)
%
% Inputs:
%    matrix - numerical matrix
%    cPZ - conPolyZonotope object 
%
% Outputs:
%    cPZ - conPolyZonotope object after multiplication with a matrix
%
% Example: 
%    c = [0;0];
%    G = [2 0 2; 0 2 2];
%    expMat = [1 0 1; 0 1 1; 0 0 0];
%    A = [2 2 4 -4];
%    b = 0;
%    expMat_ = [1 0 1 0; 0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    M = [3 1;2 -4];
%    cPZ_ = M * cPZ;
%
%    figure, hold on;
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',15);
%
%    figure; hold on;
%    plot(cPZ_,[1,2],'b','Filled',true,'EdgeColor','none','Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/mtimes, plus

% Author:       Niklas Kochdumper
% Written:      15-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % normal matrix
    if isnumeric(matrix)
        
        cPZ.c = matrix*cPZ.c;
        cPZ.G = matrix*cPZ.G;
    
        if ~isempty(cPZ.Grest)
            cPZ.Grest = matrix*cPZ.Grest;
        end
        
    % use polynomial zonotope method for interval matrices
    else
        
        pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.Grest,cPZ.expMat);
        
        temp = mtimes(matrix,pZ);
        
        cPZ = conPolyZono(temp.c,temp.G,temp.expMat,cPZ.A,cPZ.b, ...
                          cPZ.expMat_,temp.Grest,cPZ.id);
    end
end

%------------- END OF CODE --------------