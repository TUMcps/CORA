function a = set(a,varargin)
% Purpose:  Set asset properties from the specified object
% Pre:      simulation object
% Post:     property value
% Tested:   12.03.08,MA

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
    case 'initialSet'
        a.initialSet = val;        
    otherwise
        error('Asset properties: postion, speed')
    end
end