function obj = highOrderMappingMatrix(obj,intermediateOrder)
% highOrderMappingMatrix - computes a mapping matrix set without the first
% two orders
%
% Syntax:  
%    obj = highOrderMappingMatrix(obj,intermediateOrder)
%
% Inputs:
%    obj - linParamSys object 
%    intermediateOrder - order until which the original matrix set
%    representation is used
%
% Outputs:
%    obj - resulting linParamSys object 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      05-August-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%powers
zPow = obj.power.zono;
iPow = obj.power.int;

%remainder
E = obj.E;

%step size
r=obj.stepSize;

%zonotope computations
eZ = zeros(obj.dim);
eZ_input = zeros(obj.dim);
for i=3:intermediateOrder
    eZ = eZ + zPow{i}*(1/factorial(i));
    eZ_input = eZ_input + zPow{i}*(r/factorial(i+1));
end

%interval computations
eI = zeros(obj.dim);
eI_input = zeros(obj.dim);
for i=(intermediateOrder+1):obj.taylorTerms
    eI = eI + iPow{i}*(1/factorial(i));
    eI_input = eI_input + iPow{i}*(r/factorial(i+1));
end
eI = eI + E;
eI_input = eI_input + E*r;

%center of interval computations
eImid = mid(eI.int);
eImid_input = mid(eI_input.int);

%save results
obj.mappingMatrixSet.highOrderZono = eZ + eImid;
obj.mappingMatrixSet.highOrderInt = eI + (-eImid);

obj.mappingMatrixSet.highOrderZonoInput = eZ_input + eImid_input;
obj.mappingMatrixSet.highOrderIntInput = eI_input + (-eImid_input);


%------------- END OF CODE --------------