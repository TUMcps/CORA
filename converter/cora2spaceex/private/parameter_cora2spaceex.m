function parameter_cora2spaceex(Obj, component, docNode)
% parameter_cora2spaceex - 
%
% Syntax:
%    parameter_cora2spaceex(Obj, component, docNode)
%
% Inputs:
%    Obj -
%    component -
%    docNode - 
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

% get number of states and inputs
[nx,nu] = aux_numberOfParams(Obj);
params = [nx;nu];

parameters = {'x', 'u'}; % two types of parameters (x and u)

for idx_vector = 1: numel (parameters)
    for idx = 1:params(idx_vector)
            
        x = sprintf('%d', idx);
        x = [parameters{idx_vector},x];

        %Add the element node (param), for the parent element (component) and
        %set the name, type, local, d1, d2 and dynamics attribute.
        param = docNode.createElement('param');
        param.setAttribute('name',x);
        param.setAttribute('type','real');
        param.setAttribute('local', 'false');
        param.setAttribute('d1','1');
        param.setAttribute('d2','1');
        param.setAttribute('dynamics','any');
        component.appendChild(param);
    end
end


end


% Auxiliary functions -----------------------------------------------------

function [nx,nu] = aux_numberOfParams(obj)
% get the number of inputs and number of states. Only inputs that really
% act on the system dynamics and are not just dummy inputs are considered.

    if isa(obj,'nonlinearSys')
        temp = inputArgsLength(obj.mFile,2);
        nx = temp(1);
        nu = temp(2);
    elseif isa(obj,'linearSys')
        nx = size(obj.A,1);
        B = obj.B;
        B = B(:,sum(abs(B),1) > 0);
        if isscalar(B)
            nu = nx;
        else
            nu = size(B,2);
        end
    elseif isa(obj,'hybridAutomaton')
        locs = obj.location;
        nx = 0;
        nu = 0;
        for i = 1:length(locs)
           sys = locs(i).contDynamics;
           [nxTemp,nuTemp] = aux_numberOfParams(sys);
           nx = max(nx,nxTemp);
           nu = max(nu,nuTemp);
        end
    else
        throw(CORAerror('CORA:notSupported','Given contDynamics class not supported.'));
    end
end

% ------------------------------ END OF CODE ------------------------------
