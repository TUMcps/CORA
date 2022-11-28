function res = bool2stl(obj)
% bool2stl - convert boolean to stl object 
%
% Syntax:  
%    res = bool2stl(obj)
%
% Inputs:
%    obj - boolean
%
% Outputs:
%    res - resulting stl object (class stl)
%
% Example: 
%    res = bool2stl(true)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = stl('x',1);
    res.type = 'true';
end

%------------- END OF CODE --------------