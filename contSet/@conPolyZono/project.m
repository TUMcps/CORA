function res = project(obj,dim)
% projects - Projects a conPolyZonotope object onto a lower dimensional
%            subspace
%
% Syntax:  
%    res = project(obj,dim)
%
% Inputs:
%    obj - conPolyZonotope object
%    dim - dimensions spanning the subspace
%
% Outputs:
%    res - projected conPolyZonotope object
%
% Example: 
%    c = [0;0;0];
%    G = [2 0 1;0 2 1;1 0 0];
%    expMat = [1 0 3;0 1 1;0 0 0];
%    A = [1 1 -1.5];
%    b = 0.5;
%    expMat_ = [1 0 0; 0 1 0; 0 0 1];
%    Grest = [0.5 0.2; 0 0.1; 0.2 0];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest);
%     
%    cPZ_ = project(cPZ,[1,2]);
%
%    figure; hold on; box on; grid on;
%    plot(cPZ,[1,2,3],'b','Splits',10,'Filled',true);
%    plot(cPZ_,[1,2],'r','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project, conZonotope/project

% Author:       Niklas Kochdumper
% Written:      07-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = obj;
    res.c = obj.c(dim);
    res.G = obj.G(dim,:);
    
    if ~isempty(obj.Grest)
       res.Grest = res.Grest(dim,:); 
    end
end
    
%------------- END OF CODE --------------