function [obj] = tie(obj,options)
% tie - tie: time interval error; computes the error done by 
% building the convex hull of time point solutions
%
% Syntax:  
%    [obj]=tie(obj,options)
%
% Inputs:
%    obj - linearSys object
%
% Outputs:
%    obj - linearSys object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: expm, inputSol

% Author: Matthias Althoff
% Written: 08-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from object/options structure
Apower=obj.taylor.powers;
taylorTerms=options.taylorTerms;
r=options.timeStep;
dim=dimension(obj);

%initialize Asum
Asum=zeros(dim);

for i=2:taylorTerms
    %compute factor
    exp1=-i/(i-1); exp2=-1/(i-1);
    int=interval((i^exp1-i^exp2)*r^i,0);    
    %compute powers
    Aadd=int*Apower{i};
    %compute sum
    Asum=Asum+Aadd/factorial(i);
end

%write to object structure
obj.taylor.F=Asum+obj.taylor.error;

%------------- END OF CODE --------------