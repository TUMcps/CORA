function err = errOptNumIssue()
% errOptNumIssue - standardized format if a optimization solver returns
% with numerical issues
%
% Syntax:  
%    msg = errOptNumIssue()
%
% Inputs:
%    -
%
% Outputs:
%    msg - error message
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:        Victor Gassmann
% Written:       01-July-2021 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% name of operation
st = dbstack;
filename = st(2).name; % calling function name

msg = sprintf(...
    "SDP solver in %s() failed due to numerical/other issues!",...
    filename);
id = 'CORA:solverIssue';
err = MException(id,msg);
%------------- END OF CODE --------------
