function printInterval(I)
% printInterval - prints an interval such that if one executes this command
%    in the workspace, this interval would be created
%
% Syntax:
%    printInterval(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    --- 
%
% Example: 
%    I = interval([-1;-2],[2;3]);
%    printInterval(I)

% Authors:       Matthias Althoff
% Written:       12-August-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning("CORA:deprecated","function","printInterval(I)","CORA v2025.0.2","Please use printSet(I) instead.","");
printSet(I)

% ------------------------------ END OF CODE ------------------------------
