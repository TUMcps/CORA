function nlnsys = nonlinearSys(linsys)
% nonlinearSys - converts a linearSys object to a nonlinearSys object
%
% Syntax:
%    nlnsys = nonlinearSys(linsys)
%
% Inputs:
%    linsys - linearSys object
%
% Outputs:
%    nlnsys - nonlinearSys object
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
%    linsys = linearSys(A,B,c,C,D,k)
%
%    nlnsys = nonlinearSys(linsys)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-January-2023
% Last update:   02-September-2024 (MW, message for disturbance/noise matrices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% non-unit matrices for E and F currently not supported by nonlinearSys
if any(any(linsys.E)) || any(any(linsys.F))
    throw(CORAerror('CORA:notSupported', ...
        'Only all-zero disturbance/noise matrices supported for conversion.'));
end

% function handle for state equation
f = aux_funHan(linsys.A,linsys.B,linsys.c);

% check if an output equation is given
if ( (isscalar(linsys.C) && linsys.C == 1) || ...
        ( all(size(linsys.C) == [linsys.nrOfStates,linsys.nrOfStates]) ...
        && all(all(linsys.C == eye(linsys.nrOfStates))) ) ) ...
        || ~any(any(linsys.D)) || ~any(linsys.k)

    % convert output matrix to full matrix
    C = linsys.C;
    if isscalar(linsys.C) && linsys.C == 1
        C = eye(linsys.nrOfStates);
    end

    g = aux_funHan(C,linsys.D,linsys.k);

    % instantiate resulting nonlinearSys object
    nlnsys = nonlinearSys(linsys.name,f,g);

else

    % instantiate resulting nonlinearSys object
    nlnsys = nonlinearSys(linsys.name,f);

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
