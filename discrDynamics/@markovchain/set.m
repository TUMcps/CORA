function a = set(a,varargin)
% Purpose:  Set asset properties from the specified object
% Pre:      markovchain object
% Post:     property value
% Tested:   14.09.06,MA

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
    case 'field'
        a.field = val;      
    otherwise
        error('Asset properties: postion, speed')
    end
end