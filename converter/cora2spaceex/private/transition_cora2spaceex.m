function tran = transition_cora2spaceex(component, docNode, target, id, syncLabel)
% transition_cora2spaceex - 
%
% Syntax:
%    tran = transition_cora2spaceex(component, docNode, target, id, syncLabel)
%
% Inputs:
%    component -
%    docNode - 
%    target - number of target location
%    id - number of source location
%    syncLabel - synchronization label
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

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   19-June-2023 (MW, add synchronization label)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Add the element node (transition), for the parent element (component) and
%set the source and target values attribute.
tran = docNode.createElement('transition');
tran.setAttribute('source',num2str(id));
tran.setAttribute('target',num2str(target));
% only add label if non-empty
if ~isempty(syncLabel)
    tran.setAttribute('label',syncLabel);
end
component.appendChild(tran);

% ------------------------------ END OF CODE ------------------------------
