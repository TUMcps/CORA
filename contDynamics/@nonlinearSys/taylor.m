function taylor(nlnsys,varargin)
% taylor - computes symbolically the Taylor expansion of the nonlinear 
%    system; the result is stored in a m-file and passed by a handle
%
% Syntax:
%    taylor(nlnsys)
%    taylor(nlnsys,order,expPoint)
%
% Inputs:
%    nlnsys - nonlinearSys object
%    order - order of Taylor expansion (optional)
%    expPoint - expansion point (optional)
%
% Outputs:
%    -
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       06-December-2016 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain optional arguments
[order,expPoint] = setDefaultValues({6,zeros(sys.nrOfDims+sys.nrOfInputs,1)},varargin);

% create symbolic variables
[x,u] = symVariables(nlnsys,true);

% set path for reading and writing files
path = [CORAROOT 'contDynamics' filesep 'stateSpaceModels'];

% insert symbolic variables into the system equations
fdyn = nlnsys.mFile(x,u);

% compute Taylor expansion
disp('create Taylor expansion');
taylorVec = [x;u];

try
    fdyn_taylor = taylor(fdyn, taylorVec, expPoint, 'Order', order);
catch
    disp('Taylor model does not exist for this expansion point. Please try another expansion point.')
end

% write results to file
disp('create Taylor file');
createTaylorFile(fdyn_taylor,path,nlnsys.name);

% ------------------------------ END OF CODE ------------------------------
