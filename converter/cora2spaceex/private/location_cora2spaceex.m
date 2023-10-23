function location = location_cora2spaceex(Obj,component, docNode, name, idx)
% location_cora2spaceex - 
%
% Syntax:
%   location_cora2spaceex(Obj,component, docNode, name, idx)
%
% Inputs:
%    Obj -
%    component -
%    docNode - 
%    name -
%    idx -
%
% Outputs:
%    -
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isa(Obj,'hybridAutomaton')
    %Add the element node (location), for the parent element (component) and
    %set the id and name attribute.
    location = docNode.createElement('location');
    location.setAttribute('id',num2str(idx));
    location.setAttribute('name',name);
    component.appendChild(location);
else
    %Add the element node (location), for the parent element (component) and
    %set the id, name, x-y coordinates and width-height values attribute.
    location = docNode.createElement('location');
    location.setAttribute('id','1');
    location.setAttribute('name','always');
    location.setAttribute('x','320');
    location.setAttribute('y','90');
    location.setAttribute('width','300');
    location.setAttribute('height','120');
    component.appendChild(location);
end

end

% ------------------------------ END OF CODE ------------------------------
