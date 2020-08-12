function options = check_maxPolyZonoRatio(options, obj)
% check_maxPolyZonoRatio - checks if options.maxPolyZonoRatio
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    options = check_maxPolyZonoRatio(options,obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - updated options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      02-January-2020
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

    strct = 'options.polyZono';
    option = 'maxPolyZonoRatio';
    
   if isfield(options.polyZono,option)
       temp = options.polyZono.maxPolyZonoRatio;
       if ~isscalar(temp) || temp <= 0
           error(printOptionOutOfRange(obj,option,strct));
       end
   else
       options.polyZono.maxPolyZonoRatio = inf;
   end 
end

%------------- END OF CODE --------------
