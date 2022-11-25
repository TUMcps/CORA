function [obj] = tie(obj,options)
% tie - tie=time interval error; computes the error done by 
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

% Author:       Matthias Althoff
% Written:      08-May-2007 
% Last update:  10-August-2010
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object/options structure
Apower=obj.taylor.powers;
taylorTerms=options.taylorTerms;
rbyfac=options.factor;
dim=dimension(obj);

%initialize Asum
Asum_pos=zeros(dim);
Asum_neg=zeros(dim);

for i=2:taylorTerms
    %compute factor
    exp1=-i/(i-1); exp2=-1/(i-1);
    factor=(i^exp1-i^exp2)*rbyfac(i); 
    
    %init Apos, Aneg
    Apos=zeros(dim);
    Aneg=zeros(dim);
    
    %obtain positive and negative parts
    pos_ind = Apower{i}>0;
    neg_ind = Apower{i}<0;
    
    Apos(pos_ind) = Apower{i}(pos_ind);
    Aneg(neg_ind) = Apower{i}(neg_ind);
    
    %compute powers; factor is always negative
    Asum_pos=Asum_pos + factor*Aneg; 
    Asum_neg=Asum_neg + factor*Apos;
end
%instantiate interval matrix
Asum = interval(Asum_neg,Asum_pos);

%write to object structure
obj.taylor.F=Asum+obj.taylor.error;

%------------- END OF CODE --------------