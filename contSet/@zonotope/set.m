function obj = set(obj, varargin)
% set - Set data of obj
%
% Syntax:  
%    obj = set(obj, varargin)
%
% Properties:
%    Z - zonotope matrix

% Author: Matthias Althoff
% Written: 30-September-2006 
% Last update: 23-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    property = propertyArgIn{1};
    value = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch property
        case 'Z'
            obj.Z = value;        
        otherwise
            error('Property unknown')
    end
end

%------------- END OF CODE --------------