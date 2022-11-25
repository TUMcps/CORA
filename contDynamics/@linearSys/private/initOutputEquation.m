function [C,D,k] = initOutputEquation(obj,options)
% initOutputEquation - extract the C, D and k matrix
%  for the output equation y = C x + D u + k
% note: largely copied from initOutputEquation of @contDynamics > reach.m
%
% Syntax:  
%    [C,D,k] = initOutputEquation(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    [C,D,k] - matrices in y = C x + D u + k
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       26-June-2019 (from @contDynamics > reach.m)
% Last update:   ---
% Last revision: ---


%------------- BEGIN CODE --------------
   
% extract output equation parameter (y = C x + D u + k)
C = obj.C; D = obj.D; k = obj.k;
if isempty(D) && isempty(k)
    if isempty(C)
        error('All parameters for the output equation are empty!')
        % should not occur, since C initialized in linearSys
    elseif length(obj.A) ~= 1 && length(obj.C) == 1
        % C has default '1' from linearSys constructor
%         C = eye(obj.dim);
        n = obj.dim;
    else
        n = size(obj.C,1);
    end
else
    n = size(obj.D,1);
end

% initialize empty output equation parameter
if isempty(C)
    C = zeros(n,size(obj.A,1)); 
end
if isempty(D)
    Z = options.U.Z;
    D = zeros(n,size(Z,1)); 
end  
if isempty(k)
    k = zeros(n,1); 
end



end