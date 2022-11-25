function [eAt, eAtInt] = exponentialMatrices(obj,options)
% exponentialMatrices - computes the exponential matrix for the initial
% state and constant input solution
%
% Syntax:  
%    [eAt, eAtInt] = exponentialMatrices(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Outputs:
%    eAt - exponential matrix for the initial state solution
%    eAtInt - exponential matrix for the constant input solution
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      27-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% compute exponential matrix
obj = exponential(obj,options);
% compute reachable set due to input
obj = inputSolution(obj,options);

%extract data
eAt = obj.taylor.eAt;
eAtInt = obj.taylor.eAtInt;


%------------- END OF CODE --------------