function a = set(a,varargin)
% set - Set asset properties from the specified object
% Pre:      markovchain object
% Post:     property value
%
% Syntax:
%    a = set(a,varargin)
%
% Inputs:
%    ???
%
% Outputs:
%    ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       14-September-2006
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
    case 'field'
        a.field = val;      
    otherwise
        throw(CORAerror('CORA:specialError','Asset properties: postion, speed'))
    end
end

% ------------------------------ END OF CODE ------------------------------
