function obj = priv_inputSolution(obj,params,options)
% priv_inputSolution - computes the bloating due to the input 
%
% Syntax:
%    obj = priv_inputSolution(obj,params,options)
%
% Inputs:
%    obj - linParamSys object
%    params - model parameters
%    options - options struct
%
% Outputs:
%    obj - linParamSys object
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: expm, tie

% Authors:       Matthias Althoff
% Written:       06-August-2010 
% Last update:   25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%set of possible inputs
V=obj.B*params.U;

%compute vTrans if possible
try
    vTrans=obj.B*params.uTrans;
catch
    vTrans=[];
end

%load data from object/options structure
r=obj.stepSize;

%initialize the reachable set due to input
inputSet = V*r;
intM = eye(obj.nrOfDims)*r; %integral of the mapping matrix

%matrix zonotope
for i=1:length(obj.power.zono)
    taylorTerm = obj.power.zono{i}*(r/factorial(i+1));
    inputSet = inputSet + taylorTerm*V;
    intM = intM + intervalMatrix(taylorTerm);
end
%interval matrix
for i=(length(obj.power.zono)+1):length(obj.power.int)
    taylorTerm = obj.power.int{i}*(r/factorial(i+1));
    inputSet = inputSet + taylorTerm*V;
    intM = intM + taylorTerm;
end

%remainder term
Vabs=zonotope(interval(-1,1)*supremum(abs(interval(V))));
inputSet = inputSet + obj.E*r*Vabs;
intM = intM + obj.E*r;

%input solution due certain input
inputSetTrans = intM*zonotope(vTrans);

%delete zero generators in zonotope representation
inputSet=compact_(inputSet,'zeros',eps);
inputSetTrans=compact_(inputSetTrans,'zeros',eps);

%compute additional uncertainty if origin is not contained in input set
if options.originContained
    inputCorr = zeros(obj.nrOfDims,1);
else
    %compute inputF
    obj = priv_inputTie(obj,options);
    inputCorr = obj.inputF*zonotope(vTrans);
end

%write to object structure
obj.Rinput = reduce(inputSet + inputSetTrans,options.reductionTechnique,options.zonotopeOrder);
obj.Rtrans = inputSetTrans;
obj.RV = reduce(inputSet,options.reductionTechnique,options.zonotopeOrder);
obj.inputCorr = inputCorr;

% ------------------------------ END OF CODE ------------------------------
