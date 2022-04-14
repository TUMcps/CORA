function S = getSubset(cPZ,id,val)
% getSubset - extract a subset by specifying new ranges for the factors
%
% Syntax:  
%    S = getSubset(cPZ,id,val)
%
% Inputs:
%    cPZ - conPolyZono object
%    id - vector containing the identifiers of the factors that are changed
%    val - new values for the selected factors (interval or vector) 
%
% Outputs:
%    S - extracted subset (conPolyZono object or single point)
%
% Example: 
%    c = [0;0];
%    G = [2 0 2; 0 2 2];
%    expMat = [1 0 1; 0 1 1; 0 0 0];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    expMat_ = [0 1 2; 1 0 0; 0 1 0];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%    
%    S = getSubset(cPZ,2,interval(-0.5,0.3));
%
%    figure; hold on;
%    plot(S,[1,2],'b','Filled',true,'EdgeColor','none','Splits',15);
%    plot(cPZ,[1,2],'r','LineWidth',2,'Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono, polyZonotope/getSubset

% Author:       Niklas Kochdumper
% Written:      20-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % get subset for states by calling get subset method for polyZonotope
    pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.Grest,cPZ.expMat,cPZ.id);
    pZ = getSubset(pZ,id,val);

    % get subset for constraints 
    if ~isempty(cPZ.A)
        
       % construct polynomial zonotope for constraints
       temp = polyZonotope(-cPZ.b,cPZ.A,[],cPZ.expMat_,cPZ.id);
       
       % get subset for polynomial zonotope of the constraints
       conSet = getSubset(temp,id,val);
       
       % check for emptyness
       if isnumeric(conSet) 
           if ~all(abs(conSet) < 1e-8)
                S = []; return; 
           else
                S = conPolyZono(pZ); return; 
           end
       end 
       
       % add subsets for states and constraints
       if all(size(pZ.id) - size(conSet.id) == 0) && ...
          all(pZ.id -conSet.id == 0)
           S = conPolyZono(pZ.c,pZ.G,pZ.expMat,conSet.G,-conSet.c, ...
                           conSet.expMat,pZ.Grest,pZ.id);
       else
           % merge identifier vectors
           id = unique([pZ.id;conSet.id]);
           
           % adapt state exponent matrices
           expMat = pZ.expMat;
           temp = abs(length(pZ.id) - length(id));
           
           if temp > 0
               expMat = [expMat; zeros(temp,size(expMat,2))];
           end
           
           % adapt constraint exponent matrix
           expMat_ = zeros(length(id),size(conSet.expMat,2));
           
           for i = 1:length(conSet.id)
              expMat_(id == conSet.id(i),:) = conSet.expMat(i,:);
           end
           
           % construct resulting set
           S = conPolyZono(pZ.c,pZ.G,expMat,conSet.G,-conSet.c, ...
                           expMat_,pZ.Grest,id);
           
       end
       
    else
       if isnumeric(pZ)
            S = pZ;
       else
            S = conPolyZono(pZ);
       end
    end
end

%------------- END OF CODE --------------