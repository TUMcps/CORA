function text = benchmark_linear_verify_ARCH23_rand_RAND01
% benchmark_linear_verify_ARCH23_rand_RAND01 - evaluation of random
%    verification benchmark (verifiable)
%
% Syntax:
%    benchmark_linear_verify_ARCH23_rand_RAND01
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
filename = [CORAROOT filesep 'rand01.json'];
if ~exist(filename,"file")
    % fallback option: use test json file
    filename = [CORAROOT '/models/Cora/ARCH/AFF/rand01_test.json'];
end
[sys,params,spec] = json2cora_linearSys(filename);

% convert to zonotopes
params.R0 = zonotope(params.R0);
params.U = zonotope(params.U);

% our automated verification/falsification algorithm
options.verifyAlg = 'reachavoid:zonotope';

% solve verification task
timerVal = tic;
res_verify = verify(sys,params,options,spec);
tComp_verify = toc(timerVal);
disp("Verification successful? " + res_verify);

% output text
text = ['Random,RAND01,',num2str(res_verify),',',num2str(tComp_verify)];

end

% ------------------------------ END OF CODE ------------------------------
