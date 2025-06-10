function text = benchmark_linear_verify_ARCH23_rand_RAND02
% benchmark_linear_verify_ARCH23_rand_RAND02 - evaluation of random
%    verification benchmark (falsifiable)
%
% Syntax:
%    benchmark_linear_verify_ARCH23_rand_RAND02
%
% Inputs:
%    -
%
% Outputs:
%    text - char array

% Authors:       Mark Wetzlinger
% Written:       28-June-2024
% Last update:   29-June-2024 (TL, added fallback option for CI)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read from JSON file
filename = [CORAROOT filesep 'rand02.json'];
if ~exist(filename,"file")
    % fallback option: use test json file
    filename = [CORAROOT '/models/Cora/ARCH/AFF/rand02_test.json'];
end
[sys,params,spec] = json2cora_linearSys(filename);

% convert to zonotopes
params.R0 = zonotope(params.R0);
params.U = zonotope(params.U);

% our automated verification/falsification algorithm
options.verifyAlg = 'reachavoid:zonotope';

% solve verification task
timerVal = tic;
res = verify(sys,params,options,spec);
tComp_verify = toc(timerVal);
disp("Verification successful? " + res);

% output text
text = ['Random,RAND02,',num2str(res),',',num2str(tComp_verify)];

end

% ------------------------------ END OF CODE ------------------------------
