function res = test_linearSys_generateRandom
% test_linearSys_generateRandom - unit test for random generation of linear
%    systems
%
% Syntax:
%    res = test_linearSys_generateRandom
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       13-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no input arguments
linsys = linearSys.generateRandom();

% input values
n = 6;
nrInputs = 3;
nrOutputs = 2;
realInt = interval(-10,-4);
imagInt = interval(-2,3);

% state dimension given
linsys = linearSys.generateRandom('StateDimension',n);
assert(linsys.nrOfStates == n);

% input dimension given
linsys = linearSys.generateRandom('InputDimension',nrInputs);
assert(linsys.nrOfInputs == nrInputs);

% output dimension given
linsys = linearSys.generateRandom('OutputDimension',nrOutputs);
assert(linsys.nrOfOutputs == nrOutputs);

% all dimensions given
linsys = linearSys.generateRandom('StateDimension',n,...
    'InputDimension',nrInputs,'OutputDimension',nrOutputs);
assert(linsys.nrOfStates == n && linsys.nrOfInputs == nrInputs ...
    && linsys.nrOfOutputs == nrOutputs);

% real interval given
linsys = linearSys.generateRandom('realInterval',realInt);
ev = eigs(linsys.A);
assert(all(contains(realInt,real(ev)')));

% imaginary interval given
linsys = linearSys.generateRandom('ImaginaryInterval',imagInt);
ev = eigs(linsys.A);
assert(all(contains(imagInt,imag(ev)')));

% real and imaginary interval given
linsys = linearSys.generateRandom('realInterval',realInt,...
    'ImaginaryInterval',imagInt);
ev = eigs(linsys.A);
assert(all(contains(realInt,real(ev)')) && all(contains(imagInt,imag(ev)')));

% all properties given
linsys = linearSys.generateRandom('StateDimension',n,...
    'InputDimension',nrInputs,'OutputDimension',nrOutputs,...
    'realInterval',realInt,'ImaginaryInterval',imagInt);
ev = eigs(linsys.A);
assert(linsys.nrOfStates == n && linsys.nrOfInputs == nrInputs ...
    && linsys.nrOfOutputs == nrOutputs && all(contains(realInt,real(ev)')) ...
    && all(contains(imagInt,imag(ev)')));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
