function cPZ = updateConstraints(cPZ,cPZ1,cPZ2)
% updateConstraints - Update constraints after combining two constraint
%                     polynomial zonotopes
%
% Syntax:  
%    cPZ = updateConstraints(cPZ,cPZ1,cPZ2)
%
% Inputs:
%    cPZ - conPolyZono object whos constraints are updated
%    cPZ1,cPZ2 - conPolyZono objects whos constraints are combined to
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

% Author:        Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------     

     cPZ.A = blkdiag(cPZ1.A,cPZ2.A);
     cPZ.b = [cPZ1.b;cPZ2.b];

     if isempty(cPZ1.A)
         if ~isempty(cPZ2.A)
             temp = zeros(length(cPZ1.id),size(cPZ2.expMat_,2));
             cPZ.expMat_ = [temp;cPZ2.expMat_];
         end
     else
         if isempty(cPZ2.A)
             temp = zeros(length(cPZ2.id),size(cPZ1.expMat_,2));
             cPZ.expMat_ = [cPZ1.expMat_;temp];
         else
             cPZ.expMat_ = blkdiag(cPZ1.expMat_,cPZ2.expMat_);
         end
     end
end

%------------- END OF CODE --------------