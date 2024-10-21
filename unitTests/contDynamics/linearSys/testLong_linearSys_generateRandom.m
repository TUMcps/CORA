function res = testLong_linearSys_generateRandom
% testLong_linearSys_generateRandom - unit test for random
%    generation of linear systems
%
% Syntax:
%    res = testLong_linearSys_generateRandom
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

% values, ranges
states = [2;10;50];
inputs = [2;4;10];
outputs = [1;2;5];
realInts = {interval(-100,-1),interval(-1,-0.1),interval(-0.001,-0.00001)};
imagInts = {interval(-5,5),interval(-1,1),interval(-0.01,0.01)};

% number of tests per setting
nrTests = 5;

% loop over all combinations
for n=1:length(states)
for m=1:length(inputs)
for ell=1:length(outputs)
for rr=1:length(realInts)
for ii=1:length(imagInts)

    % input values
    nrStates = states(n);
    nrInputs = inputs(m);
    nrOutputs = outputs(ell);
    realInt = realInts{rr};
    imagInt = imagInts{ii};
    
    for t=1:nrTests
        seed = randi(2^32);
        rng(seed);

        % generate random linear system
        linsys = linearSys.generateRandom('StateDimension',nrStates,...
            'InputDimension',nrInputs,'OutputDimension',nrOutputs,...
            'realInterval',realInt,'ImaginaryInterval',imagInt);
        
        % compute eigenvalues
        ev = eigs(linsys.A);
        
        % check if all properties satisfy prescribed values/ranges
        assertLoop(linsys.nrOfStates == nrStates,n,m,ell,rr,ii)
        assertLoop(linsys.nrOfInputs == nrInputs,n,m,ell,rr,ii)
        assertLoop(linsys.nrOfOutputs == nrOutputs,n,m,ell,rr,ii)
        assertLoop(all(contains(realInt,real(ev)')),n,m,ell,rr,ii)
        assertLoop(all(contains(imagInt,imag(ev)')),n,m,ell,rr,ii)

    end

end
end
end
end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
