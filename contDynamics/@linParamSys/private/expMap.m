function [obj,Zres] = mappingMatrix(obj,Z,options)
% expMap - Computes the map of a zonotope with a matrix exponential whose
% matrix contains uncertain matrices
%
% Syntax:  
%    [Z] = expMap(Z,linIntSys,options)
%
% Inputs:
%    Z - zonotope object 
%    linIntSys - linear interval system
%    options - options struct
%
% Outputs:
%    Zres - resulting zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author: Matthias Althoff
% Written: 07-January-2009
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get time step, Taylor terms
r=options.timeStep;
taylorTerms=options.taylorTerms;

%get dimension of the system
A=linIntSys.A; 
dim=length(A{1});

%load data from object structure
N=linIntSys.taylor.N;
%compute middle and delta values
midN=center(N);
deltaN=N-midN;

%generate identity matrix
I=eye(dim);

%compute linear part:
%modify first uncertain matrix:
A{1}=A{1}+(I+midN)/r;

%get linear map
Zlin=r*(A*Z);

%get nonlinear map
Znonlin=deltaN*Z; %<-- middle part computed seperately anyway

%obtain overall solution
Zres=Zlin+Znonlin;

%compute modified A:
for i=1:length(A)
    modA{i}=r*A{i};
end

%write to object structure
linIntSys.taylor.deltaN=deltaN;
linIntSys.taylor.modA=modA;

%------------- END OF CODE --------------