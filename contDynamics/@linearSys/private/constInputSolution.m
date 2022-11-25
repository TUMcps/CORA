function [obj] = constInputSolution(obj,options)
% constInputSolution - computes the bloating due to constant input 
%
% Syntax:  
%    [obj] = constInputSolution(obj,options)
%
% Inputs:
%    obj - linearSys object
%
% Outputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      03-May-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


vTrans=obj.B*options.uTrans;

A=obj.A;
Apower=obj.taylor.powers;
E=obj.taylor.error;
taylorTerms=options.taylorTerms;
r=options.timeStep;
dim=length(A);


%init Asum
Asum=r*eye(dim);
%compute higher order terms
for i=1:taylorTerms   
    %compute sums
    Asum=Asum+Apower{i}*options.factor(i+1);
end

%compute solution due to constant input
eAtInt=Asum+E*r;
inputSolVtrans=eAtInt*zonotope(vTrans);

%compute additional uncertainty if origin is not contained in input set
if options.originContained
    inputCorr=zeros(dim,1);
else
    %compute inputF
    [obj]=inputTie(obj,options);
    inputF=obj.taylor.inputF;
    inputCorr=inputF*zonotope(vTrans);
end


%write to object structure
obj.taylor.Rtrans=inputSolVtrans;
obj.taylor.inputCorr=inputCorr;

%------------- END OF CODE --------------