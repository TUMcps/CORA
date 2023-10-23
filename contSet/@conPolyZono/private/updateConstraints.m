function cPZ = updateConstraints(cPZ,cPZ1,cPZ2)
% updateConstraints - Update constraints after combining two constrained
%    polynomial zonotopes
%
% Syntax:
%    cPZ = updateConstraints(cPZ,cPZ1,cPZ2)
%
% Inputs:
%    cPZ - conPolyZono object whose constraints are updated
%    cPZ1,cPZ2 - conPolyZono objects whose constraints are combined to
%                obtain the constraints of cPZ
%
% Outputs:
%    cPZ - conPolyZono object after updating the constraints
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus, quadMap

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cPZ.A = blkdiag(cPZ1.A,cPZ2.A);
cPZ.b = [cPZ1.b;cPZ2.b];

if isempty(cPZ1.A)
    if ~isempty(cPZ2.A)
        temp = zeros(length(cPZ1.id),size(cPZ2.EC,2));
        cPZ.EC = [temp;cPZ2.EC];
    end
else
    if isempty(cPZ2.A)
        temp = zeros(length(cPZ2.id),size(cPZ1.EC,2));
        cPZ.EC = [cPZ1.EC;temp];
    else
        cPZ.EC = blkdiag(cPZ1.EC,cPZ2.EC);
    end
end

% ------------------------------ END OF CODE ------------------------------
