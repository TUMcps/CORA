function sys = nonlinearSys(obj)
% nonlinearSys - converts a linearSys object to a nonlinearSys object
%
% Syntax:
%    sys = nonlinearSys(obj)
%
% Inputs:
%    obj - linearSys object
%
% Outputs:
%    sys - nonlinearSys object
%
% Example:
%    A = [-0.3780    0.2839    0.5403   -0.2962
%          0.1362    0.2742    0.5195    0.8266
%          0.0502   -0.1051   -0.6572    0.3874
%          1.0227   -0.4877    0.8342   -0.2372];
%    B = 0.25 * [-2 0 3; 2 1 0; 0 0 1; 0 -2 1];
%    c = 0.05 * [-4; 2; 3; 1];
%    C = [1 1 0 0; 0 -0.5 0.5 0];
%    D = [0 0 1; 0 0 0];
%    k = [0; 0.02];
%    obj = linearSys(A,B,c,C,D,k)
%
%    sys = nonlinearSys(obj)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% function handle for state equation
f = aux_funHan(obj.A,obj.B,obj.c);

% check if an output equation is given
if ( (isscalar(obj.C) && obj.C == 1) || ...
        ( all(size(obj.C) == [obj.dim,obj.dim]) && all(all(obj.C == eye(obj.dim))) ) ) ...
        || ~any(any(obj.D)) || ~any(obj.k)

    % convert output matrix to full matrix
    C = obj.C;
    if isscalar(obj.C) && obj.C == 1
        C = eye(obj.dim);
    end

    g = aux_funHan(C,obj.D,obj.k);

    % instantiate resulting nonlinearSys object
    sys = nonlinearSys(obj.name,f,g);

else

    % instantiate resulting nonlinearSys object
    sys = nonlinearSys(obj.name,f);

end

end


% Auxiliary functions -----------------------------------------------------

function funHan = aux_funHan(stateMatrix,inputMatrix,offset)

    % number of equations, states, and inputs
    [nrEq,nrStates] = size(stateMatrix);
    nrInputs = size(inputMatrix,2);

    % instantiate symbolic variables for states and input
    x = sym('x',[nrStates,1]);
    u = sym('u',[nrInputs,1]);

    % state equation
    eqStr = cell(nrEq,1);
    for i=1:nrEq
        % evaluate i-th row of differential equation (incl. offset)
        stateStr = string(stateMatrix(i,:) * x + offset(i));
        % replace x1 by x(1)
        for j=1:nrStates
            stateStr = strrep(stateStr,"x" + j,"x(" + j + ")");
        end
    
        inputStr = string(inputMatrix(i,:) * u);
        % replace u1 by u(1)
        for j=1:nrInputs
            inputStr = strrep(inputStr,"u" + j,"u(" + j + ")");
        end
    
        % i-th equation
        eqStr{i} = char(stateStr + " + " + inputStr);
    end

    % concatenate strings
    funHan = eval("@(x,u) [" + strjoin(eqStr,"; ") + "];");

end

% ------------------------------ END OF CODE ------------------------------
