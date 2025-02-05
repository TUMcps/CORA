function a = set(a,varargin)
% set - short description of the function
% Purpose:  Set asset properties from the specified object
% Pre:      simulation object
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
% Written:       12-March-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
    case 'initialSet'
        a.initialSet = val;        
    otherwise
        throw(CORAerror('CORA:specialError','Asset properties: postion, speed'))
    end
end

% ------------------------------ END OF CODE ------------------------------
