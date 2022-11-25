function res = quadMap(varargin)
% quadMap - computes the quadratic map of a constrained polynomial zonotope
%
% Syntax:  
%    res = quadMap(cPZ,Q)
%    res = quadMap(cPZ1,cPZ2,Q)
%
% Inputs:
%    cPZ,cPZ1,cPZ2 - conPolyZono objects
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    res - resulting set as a conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [2 1 0;2 0 1];
%    expMat = [1 1 1; 0 1 0; 0 0 1; 0 0 0];
%    A = [1 1 -0.5];
%    b = 0.5;
%    expMat_ = [0 0 0; 2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%    
%    Q{1} = [1 2; -1 2];
%    Q{2} = [-3 0; 1 1];
%
%    res = quadMap(cPZ,Q);
%
%    figure; hold on
%    plot(cPZ,[1,2],'b','Filled',true,'EdgeColor','none','Splits',20);
%
%    figure; hold on
%    plot(res,[1,2],'r','Filled',true,'EdgeColor','none','Splits',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap, zonotope/quadMap, cubMap

% Author:       Niklas Kochdumper
% Written:      19-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


    if nargin == 2                                  % quadratic map
        
        % compute quadratic map for polynomial zonotopes
        cPZ = varargin{1};
        pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.Grest,cPZ.expMat,cPZ.id);
        
        pZ = quadMap(pZ,varargin{2});
        
        % convert to constraint polynomial zonotope
        res = conPolyZono(pZ.c,pZ.G,pZ.expMat, ...
                          cPZ.A,cPZ.b,cPZ.expMat_,pZ.Grest,pZ.id);
                      
    else                                            % mixed quadrtic map
        
        % compute quadratic map for polynomial zonotopes
        cPZ1 = varargin{1}; cPZ2 = varargin{2};
        pZ1 = polyZonotope(cPZ1.c,cPZ1.G,cPZ1.Grest,cPZ1.expMat,cPZ1.id);
        pZ2 = polyZonotope(cPZ2.c,cPZ2.G,cPZ2.Grest,cPZ2.expMat,cPZ2.id);
        
        pZ = quadMap(pZ1,pZ2,varargin{3});
        
        % convert to constraint polynomial zonotope
        cPZ = conPolyZono(pZ);
        
        % update constraints
        res = updateConstraints(cPZ,cPZ1,cPZ2);
    end
end

%------------- END OF CODE --------------