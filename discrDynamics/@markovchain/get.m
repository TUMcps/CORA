function val = get(a, propName)
% Purpose:  Get asset properties from the specified object
% Pre:      markovchain object
% Post:     property value
% Tested:   14.09.06,MA
% Modified: 17.08.07,MA

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
    error([propName,' Is not a valid asset property'])
end