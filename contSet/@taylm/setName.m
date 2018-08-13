function tm = setName( tm, name )
% setName - changes the name of a variable in a Taylor model
%
% Syntax:  
%    tm = setName( tm, name )
%
% Inputs:
%    tm     - a taylm
%    name   - a string
%
% Outputs:
%    tm - a taylm
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Author:       Dmitry Grebenyuk
% Written:      31-July-2017
%
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    tm.names_of_var = {name};

end

