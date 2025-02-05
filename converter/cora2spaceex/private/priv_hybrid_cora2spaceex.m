function docNode = priv_hybrid_cora2spaceex(Obj, component, docNode)
% priv_hybrid_cora2spaceex - 
%
% Syntax:
%    docNode = priv_hybrid_cora2spaceex(Obj, component, docNode)
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

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

systems = Obj.location;
trans = [];
ids = [];

% listing first the parameters elements and then locations (with it's child
% elements [invariant and flow])
priv_parameter_cora2spaceex(Obj, component, docNode);

for idx = 1:numel(systems)
    % read out current system
    system = systems(idx);

    % read out name, invariant, flow, transitions
    name = system.name;
    invariant = system.invariant;
    contDynamics = system.contDynamics;
    trans_ = system.transition;
    
    % read out location
    location = priv_location_cora2spaceex(Obj, component, docNode, name, idx);
    
    % convert invariant and flow
    priv_invariant_cora2spaceex(Obj,location, docNode, invariant);
    priv_flow_cora2spaceex(contDynamics, location, docNode);
    
    % concatenate transitions correctly (only non-empty transitions)
    for j=1:length(trans_)
        if ~isemptyobject(trans_)
            try
                trans = [trans; trans_];
            catch
                trans = [trans; trans_'];
            end
            ids = [ids; ones(length(trans_),1)*idx];
        end
    end
    
end

% listing the transition elements with it's [guard and assignment] child
% elements
for idx = 1:numel(trans)
    % read out transition
    tran = trans(idx);

    % read out reset, guard set, target, synchronization label
    reset = tran.reset;
    guard = tran.guard;
    target = tran.target;
    syncLabel = tran.syncLabel;
    id = ids(idx);
    
    % TODO: integrate syncLabel in conversion!
    tran = priv_transition_cora2spaceex(component, docNode, target, id, syncLabel);
    priv_guard_cora2spaceex(tran,docNode,guard);    
    priv_assignment_cora2spaceex(tran, docNode, reset) 
end

% ------------------------------ END OF CODE ------------------------------
