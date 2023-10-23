function res = test_linParamSys_linParamSys
% test_linParamSys_linParamSys - unit test for linParamSys constructor
%
% Syntax:
%    res = test_linParamSys_linParamSys
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no input arguments
sys = linParamSys();
res = sys.dim == 0 && sys.nrOfInputs == 0 && sys.nrOfOutputs == 0;

% variables for system
A = [-2 0; 1.5 -3];
Aw = [0 0; 0.5 0];
A_int = intervalMatrix(A,Aw);
A_zon = matZonotope(A_int);
B = [1; 1];
Bw = [0.2; 0];
B_int = intervalMatrix(B,Bw);

% check different instantiations
sys = linParamSys(A_int,B);
res(end+1,1) = sys.dim == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && sys.constParam;
sys = linParamSys(A_zon,B);
res(end+1,1) = sys.dim == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && sys.constParam;
sys = linParamSys(A_int,B,'varParam');
res(end+1,1) = sys.dim == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && ~sys.constParam;
sys = linParamSys('system',A_int,B_int,'varParam');
res(end+1,1) = strcmp(sys.name,'system') ...
    && sys.dim == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && ~sys.constParam;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
