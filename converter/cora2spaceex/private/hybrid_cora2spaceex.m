function docNode = hybrid_cora2spaceex(Obj, component, docNode)
% hybrid_cora2spaceex - 
%
% Syntax:
%    docNode = hybrid_cora2spaceex(Obj, component, docNode)
%
% Inputs:
%    Obj -
%    component -
%    docNode - 
%
% Outputs:
%    docNode - 
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

systems = Obj.location;
trans = {}; ids = [];

% listing first the parameters elements and then locations (with it's child
% elements [invariant and flow])
parameter_cora2spaceex(Obj, component, docNode);

for idx = 1: numel (systems)
    system       = systems{idx};
    name         = system.name;
    invariant    = system.invariant;
    contDynamics = system.contDynamics;
    transition   = system.transition;
    
    location     = location_cora2spaceex(Obj, component, docNode, name, idx);
    
    % convert invariant and flow
    invariant_cora2spaceex(Obj,location, docNode, invariant);
    flow_cora2spaceex(contDynamics, location, docNode);
    
    try
        trans = [trans; transition];
    catch
        trans = [trans; transition'];
    end
    ids = [ids; ones(length(transition),1)*idx];
    
end

% listing the transition elements with it's [guard and assignment] child
% elements
for idx = 1: numel (trans)
    tran       = trans{idx} ;
    reset      = tran.reset;
    guard      = tran.guard;
    target     = tran.target;
    id         = ids(idx);
    
    tran       = transition_cora2spaceex(component, docNode, target, id);
    guard_cora2spaceex(tran,docNode,guard);    
    assignment_cora2spaceex(tran, docNode, reset)
    
end

end

%------------- END OF CODE --------------