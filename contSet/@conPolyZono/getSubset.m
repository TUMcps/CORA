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
%    E = [1 0 1; 0 1 1; 0 0 0];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    EC = [0 1 2; 1 0 0; 0 1 0];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    
%    S = getSubset(cPZ,2,interval(-0.5,0.3));
%
%    figure; hold on;
%    plot(S,[1,2],'FaceColor','b','Splits',12);
%    plot(cPZ,[1,2],'r','LineWidth',2,'Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono, polyZonotope/getSubset

% Authors:       Niklas Kochdumper
% Written:       20-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{cPZ,'att','conPolyZono'};
                {id,'att','numeric',{'vector','integer','nonnegative'}};
                {val,'att',{'numeric','interval'},'vector'}});

% get subset for states by calling get subset method for polyZonotope
pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E,cPZ.id);
pZ = getSubset(pZ,id,val);

% get subset for constraints 
if ~isempty(cPZ.A)
    
   % construct polynomial zonotope for constraints
   conPZ = polyZonotope(-cPZ.b,cPZ.A,[],cPZ.EC,cPZ.id);
   
   % get subset for polynomial zonotope of the constraints
   subConPZ = getSubset(conPZ,id,val);
   
   % check for emptyness
   if isnumeric(subConPZ) 
       if ~all(abs(subConPZ) < 1e-8)
            S = []; return; 
       else
            S = conPolyZono(pZ); return; 
       end
   end 
   
   % add subsets for states and constraints
   if all(size(pZ.id) - size(subConPZ.id) == 0) && ...
      all(pZ.id -subConPZ.id == 0)
       S = conPolyZono(pZ.c,pZ.G,pZ.E,subConPZ.G,-subConPZ.c, ...
                       subConPZ.E,pZ.GI,pZ.id);
   else
       % merge identifier vectors
       id = unique([pZ.id;subConPZ.id]);
       
       % adapt state exponent matrices
       E = pZ.E;
       conPZ = abs(length(pZ.id) - length(id));
       
       if conPZ > 0
           E = [E; zeros(conPZ,size(E,2))];
       end
       
       % adapt constraint exponent matrix
       E_ = zeros(length(id),size(subConPZ.E,2));
       
       for i = 1:length(subConPZ.id)
          E_(id == subConPZ.id(i),:) = subConPZ.E(i,:);
       end
       
       % construct resulting set
       S = conPolyZono(pZ.c,pZ.G,E,subConPZ.G,-subConPZ.c,E_,pZ.GI,id);
       
   end
   
else
   if isnumeric(pZ)
        S = pZ;
   else
        S = conPolyZono(pZ);
   end
end

% ------------------------------ END OF CODE ------------------------------
