function obj = inputSolution(obj,options)
% inputSolution - computes the bloating due to the input 
%
% Syntax:
%    obj = inputSolution(obj,options)
%
% Inputs:
%    obj - linProbSys object
%    options - options struct
%
% Outputs:
%    obj - linProbSys object
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
% Written:       13-February-2007 
% Last update:   26-February-2008
%                08-September-2009
%                17-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%set of possible inputs
V = obj.B*options.U;

%compute vTrans
vTrans = obj.B*options.uTrans;

%load data from object/options structure
A = obj.A;
Apower = obj.taylor.powers;
E = obj.taylor.error;
r = options.timeStep;
n = obj.dim;
I = eye(n);
F = obj.taylor.F; 

%non-probabilistic solution--------------------------------
%init Vsum
Vsum=r*V;
%compute higher order terms
for i=1:options.taylorTerms
    %compute factor
    factor=1/factorial(i+1);    

    %compute sums
    Vsum=Vsum+Apower{i}*factor*r^(i+1)*V;
end

%compute overall solution
inputSolV=Vsum+E*r*V;

%compute solution due to constant input
inputSolVtrans=inv(A)*(expm(A*r)-I)*vTrans; % ... zonotope(vTrans) ?

%compute additional uncertainty if origin is not contained in input set
if options.originContained
    inputCorr=zeros(n,1);
else
    inputCorr=inv(A)*F*zonotope(vTrans);
end

%write to object structure
obj.taylor.V=V;
uncertainMean=inputSolV+inputSolVtrans;
obj.taylor.Rinput=probZonotope([uncertainMean.c,uncertainMean.G],zeros(n,1),options.gamma);
obj.taylor.Rtrans=inputSolVtrans;
obj.taylor.inputCorr=inputCorr;
%----------------------------------------------------------

%probabilistic solution------------------------------------
%obtain covariance matrix after one time step
[V,W]=eig(A);
lambda=diag(W); %eigenvalues

C=inv(V)*obj.C;
D=C*C.';

%compute Sigma in transformed coordinates
for i=1:length(lambda)
    for j=1:length(lambda)
        lambdaSum=lambda(i)+lambda(j);
        Sigma(i,j)=D(i,j)/lambdaSum*(1-exp(-lambdaSum*r));
    end
end
%Sigma in original coordinates
Sigma=V*Sigma*V.';
G=aux_generators(Sigma);

%instantiate probabilistic zonotope
ProbInputSol=probZonotope(zeros(n,1),G,options.gamma);

%write to object structure
obj.taylor.pRinput=ProbInputSol;

end


% Auxiliary functions -----------------------------------------------------

function [G]=aux_generators(Sigma)
% returns the generator matrix of a probabilistic zonotope if
% the covariance matrix Sigma is given

    %ensure symmetry for numerical stability
    Sigma=0.5*(Sigma+Sigma');
    
    %get eigenvalue, eigenvectors of Sigma
    [V,W]=eig(Sigma);
    
    %compute new generators
    G=V*sqrt(W);

end

% ------------------------------ END OF CODE ------------------------------
