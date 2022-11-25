function [Rfirst] = coreReach(obj,Rinit,options)
% coreReach - computes the reachable continuous set for the next time step
% without considering uncertain inputs
%
% Syntax:  
%    [Rfirst] = coreReach(obj,Rinit,options)
%
% Inputs:
%    Rinit - initial reachable set
%
% Outputs:
%    obj - linearSys object
%    Rfirst - first reachable set 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-August-2011 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%first time step homogeneous solution
Rfirst = obj.mappingMatrixSet.zono*Rinit + obj.mappingMatrixSet.int*Rinit + obj.Rtrans;

% %load matrix zonotope
% matZ = obj.mappingMatrixSet.zono;
% 
% %obtain order vector for reduction
% nEval(1) = norm(matZ.center);
% for i=1:matZ.gens
%     nEval(i+1) = norm(matZ.generator{i});
% end
% 
% %normalize nEval
% nEval = nEval/sum(nEval);
% 
% %obtain order vector
% orderVec = options.zonotopeOrder*nEval;
% for i=1:length(orderVec)
%     if orderVec(i) < 1
%         orderVec(i) = 1;
%     end
% end
% 
% Rred = reduce(Rinit,'girardMultiple',orderVec);
% 
% %zonotope matrix computation
% Rzono = matZ.center*Rred{1};
% for i = 1:matZ.gens
%     Rzono = Rzono + matZ.generator{i}*Rred{i+1};
% end
% 
% Rfirst = Rzono + obj.mappingMatrixSet.int*Rred{end} + obj.Rtrans;


%------------- END OF CODE --------------