function taylor(obj,varargin)
% taylor - computes symbolically the Taylor expansion of the nonlinear 
% system; the result is stored in a m-file and passed by a handle
%
% Syntax:  
%    [obj] = taylor(obj)
%
% Inputs:
%    obj - nonlinear system object
%    order - order of Taylor expansion (optional)
%    expPoint - expansion point (optional)
%
% Outputs:
%    ---
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      06-December-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create symbolic variables
[x,u] = symVariables(obj,'LRbrackets');

% obtain optional arguments
if nargin == 2
    order = varargin{1};
    expPoint = zeros(length([x;u]),1);
elseif nargin == 3
    order = varargin{1};
    expPoint = varargin{2};
else
    order = 6;
    expPoint = zeros(length([x;u]),1);
end

% set path for reading and writing files
path = [coraroot '/contDynamics/stateSpaceModels'];

% insert symbolic variables into the system equations
fdyn = obj.mFile(x,u);

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
createTaylorFile(fdyn_taylor,path,obj.name);

%------------- END OF CODE --------------