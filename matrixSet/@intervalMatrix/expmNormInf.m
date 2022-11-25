function val = expmNormInf(intMat,t)
% expmNormInf - returns the over-approximation of the inf-norm of the 
% symmetric solution of the computation of the mapping mimicing the 
% exponential
%
% Syntax:  
%    val = expmNormInf(intMat,t)
%
% Inputs:
%    intMat - interval matrix
%    t - time
%
% Outputs:
%    val - value of the norm
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      25-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%extract nominal matrix and symmetric interval matrix
infA = infimum(intMat.int);
supA = supremum(intMat.int);
A = 0.5*(supA+infA);
S = 0.5*(supA-infA);

%set starting order
j=2;

%derive norm values
nA = norm(A,inf);
nS = norm(S,inf);
nAS = norm(abs(A) + S,inf);

%compute summand1
summand1 = (nA*nS*t^2/j)/(nA*(nAS)*t^2/j^2 - (nA+nAS)*t/j + 1);
summand2 = nS*t/(1 - nAS*t/j);

%final result
val = summand1 + summand2;

% %Tests-----------------------------------
% %build G matrix
% G = [nA*t/j 0 0;...
%      nS*t/j nAS*t/j 0;...
%      0 1 1];
%  
%  %compute approximated Ginf
%  Ginf = G^100;
%  
%  %set initial vector
%  initVec = [nA*t; nS*t; 0];
%  
%  %final result
%  val2 = [0 0 1]*Ginf*initVec;

%------------- END OF CODE --------------