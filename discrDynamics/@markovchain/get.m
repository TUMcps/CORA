function val = get(a, propName)
% get - get asset properties from the specified object
% Pre:      markovchain object
% Post:     property value
%
% Syntax:
%    val = get(a, propName)
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
% Last update:   17-August-2007
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

switch propName
    case 'field'
        val = a.field;     
    case 'T'
        val = a.T;            
    case 'nrOfSegments'
        field = a.field;
        val = get(field,'nrOfSegments');    
    case 'actualSegmentNr'
        field = a.field;
        val = get(field,'actualSegmentNr');            
otherwise
    throw(CORAerror('CORA:specialError',[propName,' Is not a valid asset property']))
end

% ------------------------------ END OF CODE ------------------------------
