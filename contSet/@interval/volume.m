function vol = volume(I)
% volume - Computes volume of an interval
%
% Syntax:  
%    vol = volume(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    vol - volume
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    vol = volume(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      24-July-2016 
% Last update:  18-August-2022 (MW, include standardized preprocessing)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_volume('interval',I);

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    vol = vars{1}; return
end

% simple volume formula
vol = prod(2*rad(I));

%------------- END OF CODE --------------