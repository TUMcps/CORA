function val = get(obj, propName)
% get - Retrieve data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    propName - name of the property ('simulation', 'reachableSet' or 
%               'locationReach')
%
% Outputs:
%    val - value of the property
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schöpfer, Niklas Kochdumper
% Written:      08-June-2018  
% Last update:  09-July-2018
% Last revision: ---

%------------- BEGIN CODE --------------

switch propName
    case 'reachableSet'
        val = obj.result.reachSet;
    case 'simulation'
        val = obj.result.simulation; 
    case 'locationReach'
        val = cell(length(obj.result.reachSet),1);
        for i = 1:length(val)
           val{i} = obj.result.reachSet{i}.location; 
        end 
otherwise
    error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------
