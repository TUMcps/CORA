function vol = volume(P)
% volume - Computes the volume of a mptPolytope
%
% Syntax:  
%    vol = volume(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    vol - volume
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      02-February-2011 
% Last update:  18-August-2022 (MW, include standardized preprocessing)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_volume('mptPolytope',P);

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    vol = vars{1}; return
end

%call volume operation of mpt toolbox
vol = volume(P.P);

%------------- END OF CODE --------------