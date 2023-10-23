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
% See also: conZonotope/compact_
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       04-January-2020
% Last update:   28-March-2022 (MW, fix bugs for special cases)
%                22-May-2022
% Last revision: 29-July-2023 (MW, merged to compact)

% ------------------------------ BEGIN CODE -------------------------------

funcname = mfilename;
warning(sprintf(['The function ''' funcname ''' is deprecated (since CORA 2024) and has been replaced by ''compact''.\n' ...
    '         When updating the code, please rename every function call ''' funcname '(cZ)'' -> ''compact(cZ,''zeros'')''.\n' ...
    '         Note that the function ''' funcname ''' will be removed in a future release.']));
cZ = compact_(cZ,'zeros',eps);

% ------------------------------ END OF CODE ------------------------------
