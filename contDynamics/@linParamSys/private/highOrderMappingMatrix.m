function obj = highOrderMappingMatrix(obj,intermediateTerms)
% highOrderMappingMatrix - computes a mapping matrix set without the first
%    two orders
%
% Syntax:
%    obj = highOrderMappingMatrix(obj,intermediateTerms)
%
% Inputs:
%    obj - linParamSys object 
%    intermediateTerms - order until which the original matrix set
%                        representation is used
%
% Outputs:
%    obj - resulting linParamSys object 
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       05-August-2010
% Last update:   15-February-2021 (MW, rename: intermediateTerms)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
for i=3:intermediateTerms
    eZ = eZ + zPow{i}*(1/factorial(i));
    eZ_input = eZ_input + zPow{i}*(r/factorial(i+1));
end

%interval computations
eI = zeros(obj.dim);
eI_input = zeros(obj.dim);
for i=(intermediateTerms+1):obj.taylorTerms
    eI = eI + iPow{i}*(1/factorial(i));
    eI_input = eI_input + iPow{i}*(r/factorial(i+1));
end
eI = eI + E;
eI_input = eI_input + E*r;

%center of interval computations
eImid = center(eI.int);
eImid_input = center(eI_input.int);

%save results
obj.mappingMatrixSet.highOrderZono = eZ + eImid;
obj.mappingMatrixSet.highOrderInt = eI + (-eImid);

obj.mappingMatrixSet.highOrderZonoInput = eZ_input + eImid_input;
obj.mappingMatrixSet.highOrderIntInput = eI_input + (-eImid_input);

% ------------------------------ END OF CODE ------------------------------
