function obj = bool2stl(bool)
% bool2stl - convert boolean to stl object 
%
% Syntax:  
%    obj = stl.bool2stl(bool)
%
% Inputs:
%    bool - true/false
%
% Outputs:
%    obj - resulting stl object (class stl)
%
% Example: 
%    obj = stl.bool2stl(true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      09-November-2022 
% Last update:  13-April-2023 (MW, make static function, update code)
% Last revision:---

%------------- BEGIN CODE --------------

% ensure that bool is really a logical
inputArgsCheck({{bool,'att','logical','scalar'}});

obj = stl('x',1);
if bool
    obj.type = 'true';
else
    obj.type = 'false';
end
obj.variables = [];
obj.logic = true;

%------------- END OF CODE --------------