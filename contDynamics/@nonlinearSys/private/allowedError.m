function [allErr] = allowedError(obj,linSys,options)
% allowedError - computes the allowed linearization error
%
% Syntax:  
%    [allErr] = allowedError(obj,linSys,options)
%
% Inputs:
%    obj - nonlinear system object
%    linSys - linear system object
%    options - options struct
%
% Outputs:
%    obj - nonlinear system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:        Matthias Althoff
% Written:       29-October-2007 
% Last update:   22-January-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%retrieve data
f0=obj.linError.f0;
deltaT=options.timeStep;
expFactor=obj.expFactor;

A=linSys.A;
eAt=expm(A*deltaT);
dim=length(f0);
I=eye(dim);

%compute allowed lin error
%allErr=inv(eAt-I)*A*f0*deltaT*expFactor;

%compute inverse and check if possible
if det(eAt-I)==0
    inverse=inv(eAt-I+1e-4);
    disp('inversion approximation');
else
    inverse=inv(eAt-I);
end
%allErr=inverse*A*ones(dim,1)*deltaT*expFactor;
allErr=inverse*A*deltaT*expFactor;
allErr=abs(allErr);

%------------- END OF CODE --------------