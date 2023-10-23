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
% See also: ---

% Authors:       Matthias Althoff
% Written:       24-January-2007 
% Last update:   22-March-2007
%              1 09-January-2009 (vnorm acceleration)
%              1 01-October-2017 (use of auxiliary function pickedGenerators)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize Z_red
Zred=Z;

% pick generators to reduce
[center, Gunred, Gred] = pickedGeneratorsFast(Z,order);

% box remaining generators
d=sum(abs(Gred),2);
Gbox=diag(d);

%build reduced zonotope
Zred.c = center;
Zred.G = [Gunred,Gbox];


% ------------------------------ END OF CODE ------------------------------
