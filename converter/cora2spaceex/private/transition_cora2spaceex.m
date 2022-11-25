function tran = transition_cora2spaceex(component, docNode, target, id)
% transition_cora2spaceex - 
%
% Syntax:
%    tran = transition_cora2spaceex(component, docNode, target, id)
%
% Inputs:
%    component -
%    docNode - 
%    target -
%    id -
%
% Outputs:
%    tran - 
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author:        Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---
%
%------------- BEGIN CODE --------------

%Add the element node (transition), for the parent element (component) and
%set the source and target values attribute.
tran = docNode.createElement('transition');
tran.setAttribute('source',num2str(id));
tran.setAttribute('target',num2str(target));
component.appendChild(tran);

end

%------------- END OF CODE --------------