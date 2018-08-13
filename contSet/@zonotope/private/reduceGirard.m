function Zred = reduceGirard(Z,order)
% reduceGirard - Reduce zonotope so that its order stays below a specified
% limit 
%
% Syntax:  
%    [Zred]=reduceGirard(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Matthias Althoff
% Written: 24-January-2007 
% Last update: 22-March-2007
%              19-January-2009 (vnorm acceleration)
%              11-October-2017 (use of auxiliary function pickedGenerators)
% Last revision: ---

%------------- BEGIN CODE --------------

% initialize Z_red
Zred=Z;

% pick generators to reduce
try
    [center, Gunred, Gred] = pickedGenerators(Z,order);
catch
    disp('stop')
end

% box remaining generators
d=sum(abs(Gred),2);
%build box Gbox from interval hull vector d
Gbox=diag(d);

%build reduced zonotope
Zred.Z=[center,Gunred,Gbox];


%------------- END OF CODE --------------