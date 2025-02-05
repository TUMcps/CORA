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
assert(sys.nrOfDims == 0 && sys.nrOfInputs == 0 && sys.nrOfOutputs == 0);

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
assert(sys.nrOfDims == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && sys.constParam);
sys = linParamSys(A_zon,B);
assert(sys.nrOfDims == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && sys.constParam);
sys = linParamSys(A_int,B,'varParam');
assert(sys.nrOfDims == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && ~sys.constParam);
sys = linParamSys('system',A_int,B_int,'varParam');
assert(strcmp(sys.name,'system') ...
    && sys.nrOfDims == 2 && sys.nrOfInputs == 1 ...
    && sys.nrOfOutputs == 0 && ~sys.constParam);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
