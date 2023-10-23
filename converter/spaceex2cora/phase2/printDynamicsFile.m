function printDynamicsFile(path,name,dynamics,usecase)
% printDynamicsFile - create a file defining a nonlinear function for the
%    spaceex2cora converter 
%
% Syntax:
%    printDynamicsFile(path,name,dynamics,usecase)
%
% Inputs:
%    path - path to the directory created during conversion of the model
%    name - name of the generated .m file (specified as string)
%    dynamics - symbolic equations defining the function to be created
%               check eq2linSys for specifications on the format of these
%               symbolic equations
%    usecase - string specifiying the use of the created function
%              available options: flow, reset
%   
% Outputs:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       ???, Maximilian Perschl
% Written:       ---
% Last update:   30-January-2022 (MP, nonlinear reset functions in sx2cora)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~contains(['flow','reset'],usecase)
    throw(CORAerror('CORA:converterIssue','Unknown usecase for dynamics file.'));
end

auxpath = [path '/auxiliary'];
if ~exist(auxpath,'dir')
    mkdir(auxpath);
end

addpath(auxpath);

fname = strcat(auxpath,'/',name,'.m');
file = fopen(fname,'w');

if file<0
    throw(CORAerror('CORA:converterIssue',['Could not open file ' fname '.']));
end

% write file contents
if strcmp(usecase,"flow")
    Str = "function [dx]=" + name + "(x,u)" + newline + newline + dynamics;
else
    % rewrite dynamics from dx = ... to x_R = ...
    dynamics = replace(dynamics,"dx","x_R");
    % During transition.compDerivatives (called when a transition object is
    % initialized), the dynamic function is called with symbolic input.
    % If the value of the first dimension of the dynamics is hardcoded to
    % be 0, the return value of the dynamic function will be assigned
    % double which then leads to an error with symbolic inputs.
    % To circumvent this, we check the input type and assign 0 as symbolic
    % for symbolic inputs and as double for double inputs
    type_safety_str = "% Assign type of input to output variable, type safety feature" ...
                        + newline + "if isa(x,'double')" + newline ...
                        + sprintf("\tx_R(1,1) = 0;") + newline + "else" + newline + sprintf("\tx_R(1,1) = sym(0);") ...
                        + newline + "end" + newline;
                    
    Str = "function [x_R]=" + name + "(x,u)" + newline + newline + type_safety_str ...
          + newline +  newline + "% reset assignments:" + newline + dynamics;
end
fwrite(file,Str);
fclose(file);

% ensure matlab detects new function
rehash path;

% ------------------------------ END OF CODE ------------------------------
