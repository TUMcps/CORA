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
res = true;

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
if linSys.dim ~= n
    res = false;
end

% input dimension given
linSys = linearSys.generateRandom('InputDimension',nrInputs);
if linSys.nrOfInputs ~= nrInputs
    res = false;
end

% output dimension given
linSys = linearSys.generateRandom('OutputDimension',nrOutputs);
if linSys.nrOfOutputs ~= nrOutputs
    res = false;
end

% all dimensions given
linSys = linearSys.generateRandom('StateDimension',n,...
    'InputDimension',nrInputs,'OutputDimension',nrOutputs);
if linSys.dim ~= n || linSys.nrOfInputs ~= nrInputs || ...
        linSys.nrOfOutputs ~= nrOutputs
    res = false;
end

% real interval given
linSys = linearSys.generateRandom('realInterval',realInt);
ev = eigs(linSys.A);
if ~all(contains(realInt,real(ev)'))
    res = false;
end

% imaginary interval given
linSys = linearSys.generateRandom('ImaginaryInterval',imagInt);
ev = eigs(linSys.A);
if ~all(contains(imagInt,imag(ev)'))
    res = false;
end

% real and imaginary interval given
linSys = linearSys.generateRandom('realInterval',realInt,...
    'ImaginaryInterval',imagInt);
ev = eigs(linSys.A);
if ~all(contains(realInt,real(ev)')) || ~all(contains(imagInt,imag(ev)'))
    res = false;
end

% all properties given
linSys = linearSys.generateRandom('StateDimension',n,...
    'InputDimension',nrInputs,'OutputDimension',nrOutputs,...
    'realInterval',realInt,'ImaginaryInterval',imagInt);
ev = eigs(linSys.A);
if linSys.dim ~= n || linSys.nrOfInputs ~= nrInputs || ...
        linSys.nrOfOutputs ~= nrOutputs || ...
        ~all(contains(realInt,real(ev)')) || ~all(contains(imagInt,imag(ev)'))
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
