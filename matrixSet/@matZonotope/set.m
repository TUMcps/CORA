function obj = set(obj, varargin)
% set - Set data of obj
%
% Syntax:  
%    obj = set(obj, varargin)
%
% Properties:
%    matZ - matrix zonotope

% Author:       Matthias Althoff
% Written:      05-August-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    property = propertyArgIn{1};
    value = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch property
        case 'power'
            obj.power = value;        
        otherwise
            error('Property unknown')
    end
end

%------------- END OF CODE --------------