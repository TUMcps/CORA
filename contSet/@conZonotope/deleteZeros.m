function cZ = deleteZeros(cZ)
% deleteZeros - delete all constraints of the form 0 * beta = 0 as they are
%    trivially true for all values of beta; additionally, remove all
%    generators which have zero-length and no corresponding entries in the
%    constraint matrix
%    note that there are also other constraints which might be true for all
%    values of beta (-1 to 1), but this would require much more
%    computational effort and is thus omitted in this function
%
% Syntax:  
%    cZ = deleteZeros(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1 0 0 1;0 1 0 2 -1];
%    A = [-2 0 1 -1; 0 0 0 0]; b = [2;0];
%    cZ = conZonotope(Z,A,b);
%    
%    cZ_ = deleteZeros(cZ);
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(cZ_,[1,2],'b--');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/deleteZeros
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      04-January-2020
% Last update:  28-March-2022 (MW, fix bugs for special cases)
%               22-May-2022
% Last revision:---

%------------- BEGIN CODE --------------

% constraints that are trivially true: 0 * beta = 0
idx = ~any([cZ.A,cZ.b],2);
cZ.A(idx,:) = []; cZ.b(idx) = [];

% zero-length generators (corresponding columns in constraint matrix A need
% to be all-zero as well)
idx = ~any([cZ.Z(:,2:end);cZ.A],1);
% remove zero-length generators and corresponding constraints
cZ.Z(:,[false,idx]) = [];
if ~isempty(cZ.A)
    cZ.A(:,idx) = [];
end

%------------- END OF CODE --------------