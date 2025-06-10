function [linsys,params,options] = priv_initReach_Krylov(linsys,params,options)
% priv_initReach_Krylov - computes necessary preliminary results for
% reachability analysis in the Krylov subspace
%
% Syntax:
%    [linsys,params,options] = priv_initReach_Krylov(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Maximilian Perschl
% Written:       15-November-2016
% Last update:   24-July-2020 (box inputs removed)
%                25-April-2025 (MP, major refactor including new error bound for subspaces)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% augment options
options.tFinal = params.tFinal;

% initialize field for tp solutions
linsys.krylov.Rhom_tp_prev = linsys.C*params.R0;


% INPUT SOLUTION-----------------------------------------------------------

% project input set and vector
params.U = linsys.B*params.U;
uTrans = linsys.B*params.uTrans;
% compute solution for constant inputs
linsys = priv_krylov_R_uTrans(linsys,options.timeStep,uTrans,options);
% Computation of first input solution of input set
% (is necessary for subspaces)

linsys = priv_inputSolution_Krylov(linsys,params,options);

% STATE SOLUTION-----------------------------------------------------------

% compute state subspaces and save them in obj.krylov.state
linsys = aux_create_state_subspaces(linsys,params.R0,options);

% compute and save subspaces for input solution
linsys = aux_create_input_subspaces(linsys,options);

end
%------------- END OF FUNCTION PRIV_INIT_REACH_KRYLOV ---------------------


% Auxiliary functions -----------------------------------------------------

function linsys = aux_create_state_subspaces(linsys,R0,options)
% Create Krylov subspaces for initial set Rinit and save them in obj.krylov.state

% compute subspaces
[V_c,H_c,V_g,H_g] = priv_subspace_Krylov_jaweckiBound(linsys.A,R0,options);
% [V_c,H_c,V_g,H_g] = subspace_Krylov(linsys,params.R0,params,options);

% save subspace matrices within linearSys object to create taylor object and project them
% center
C = linsys.C;

if ~isempty(V_c)
    linsys.krylov.state.c_sys = linearSys(H_c,V_c');
    linsys.krylov.state.c_sys.taylor = taylorLinSys(H_c);
    linsys.krylov.state.c_sys_proj = linearSys(H_c,(C*V_c)');
    linsys.krylov.state.c_sys_proj.taylor = taylorLinSys(H_c);
else
    linsys.krylov.state.c_sys = [];
    linsys.krylov.state.c_sys_proj = [];
end

%generators
if ~isempty(V_g)
    for i=1:length(V_g)
        linsys.krylov.state.g_sys{i} = linearSys(H_g{i},V_g{i}');
        linsys.krylov.state.g_sys{i}.taylor = taylorLinSys(H_g{i});
        linsys.krylov.state.g_sys_proj{i} = linearSys(H_g{i},(C*V_g{i})');
        linsys.krylov.state.g_sys_proj{i}.taylor = taylorLinSys(H_g{i});
    end
else
    linsys.krylov.state.g_sys = [];
    linsys.krylov.state.g_sys_proj = [];
end

end
%------------- END OF FUNCTION AUX_CREATE_STATE_SUBSPACES --------------

function linsys = aux_create_input_subspaces(linsys,options)
% Create Krylov subspaces for input set U_0
% and save them in linsys.krylov.input


%% U_0
% compute subspaces
[V_c,H_c,V_g,H_g] = priv_subspace_Krylov_jaweckiBound(linsys.A,linsys.krylov.RV_0,options);
% [V_c,H_c,V_g,H_g] = subspace_Krylov(linsys,linsys.krylov.RV_0,params,options);

% save subspace matrices within linearSys object to create taylor object and project them
% center
C = linsys.C;

if ~isempty(V_c)
    linsys.krylov.input.c_sys_proj = linearSys(H_c,(C*V_c)');
    linsys.krylov.input.c_sys_proj.taylor = taylorLinSys(H_c);
else
    linsys.krylov.input.c_sys_proj = [];
end

%generators
if ~isempty(V_g)
    for i=1:length(V_g)
        linsys.krylov.input.g_sys_proj{i} = linearSys(H_g{i},(C*V_g{i})');
        linsys.krylov.input.g_sys_proj{i}.taylor = taylorLinSys(H_g{i});
    end
else
    linsys.krylov.input.g_sys_proj = [];
end

end
%------------- END OF FUNCTION AUX_CREATE_INPUT_SUBSPACES -----------------

% ------------------------------ END OF CODE ------------------------------
