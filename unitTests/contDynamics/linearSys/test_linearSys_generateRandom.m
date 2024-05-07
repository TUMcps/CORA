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

% init result
res = true(0);

% no input arguments
linSys = linearSys.generateRandom();

% input values
n = 6;
nrInputs = 3;
nrOutputs = 2;
realInt = interval(-10,-4);
imagInt = interval(-2,3);

% state dimension given
linSys = linearSys.generateRandom('StateDimension',n);
res(end+1,1) = linSys.dim == n;

% input dimension given
linSys = linearSys.generateRandom('InputDimension',nrInputs);
res(end+1,1) = linSys.nrOfInputs == nrInputs;

% output dimension given
linSys = linearSys.generateRandom('OutputDimension',nrOutputs);
res(end+1,1) = linSys.nrOfOutputs == nrOutputs;

% all dimensions given
linSys = linearSys.generateRandom('StateDimension',n,...
    'InputDimension',nrInputs,'OutputDimension',nrOutputs);
res(end+1,1) = linSys.dim == n && linSys.nrOfInputs == nrInputs ...
    && linSys.nrOfOutputs == nrOutputs;

% real interval given
linSys = linearSys.generateRandom('realInterval',realInt);
ev = eigs(linSys.A);
res(end+1,1) = all(contains(realInt,real(ev)'));

% imaginary interval given
linSys = linearSys.generateRandom('ImaginaryInterval',imagInt);
ev = eigs(linSys.A);
res(end+1,1) = all(contains(imagInt,imag(ev)'));

% real and imaginary interval given
linSys = linearSys.generateRandom('realInterval',realInt,...
    'ImaginaryInterval',imagInt);
ev = eigs(linSys.A);
res(end+1,1) = all(contains(realInt,real(ev)')) && all(contains(imagInt,imag(ev)'));

% all properties given
linSys = linearSys.generateRandom('StateDimension',n,...
    'InputDimension',nrInputs,'OutputDimension',nrOutputs,...
    'realInterval',realInt,'ImaginaryInterval',imagInt);
ev = eigs(linSys.A);
res(end+1,1) = linSys.dim == n && linSys.nrOfInputs == nrInputs ...
    && linSys.nrOfOutputs == nrOutputs && all(contains(realInt,real(ev)')) ...
    && all(contains(imagInt,imag(ev)'));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
