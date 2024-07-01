function text = benchmark_linear_verify_ARCH24_rand_generation
% benchmark_linear_verify_ARCH24_rand_generation - here, we generate the
%   random verification benchmarks for the ARCH competition
%
% Syntax:
%    benchmark_linear_verify_ARCH24_rand_generation
%
% Inputs:
%    -
%
% Outputs:
%    text - char array

% Authors:       Mark Wetzlinger
% Written:       07-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% fix seed
rng(128256512);

% dimensions
n = 10;
m = 2;
r = 3;
realInt = interval(-5,-1);
imagInt = interval(-0.5,0.5);

% specifications
nrSpecs = 2;
setRepInitialSet = 'interval';
setRepInputSet = 'interval';
setRepSpec = 'interval';
relerror = 0.1;

% our automated verification/falsification algorithm
options.verifyAlg = 'reachavoid:zonotope';


% 1. verification
sat = true;

[sys_verify,params_verify,spec_verify,sat] = generateRandomVerify(...
    "StateDimension",n,...
    "InputDimension",m,...
    "OutputDimension",r,...
    "RealInterval",realInt,...
    "ImaginaryInterval",imagInt,...
    "NrSpecs",nrSpecs,...
    "SetRepInitialSet",setRepInitialSet,...
    "SetRepInputSet",setRepInputSet,...
    "SetRepSpec",setRepSpec,...
    "Satisfiable",sat, ...
    "RelError",relerror);

% write benchmark to JSON format
benchmark2json(sys_verify,params_verify,spec_verify,'rand01.json');

% convert to zonotopes for verification algorithm
params_verify.R0 = zonotope(params_verify.R0);
params_verify.U = zonotope(params_verify.U);

% solve verification task
tic;
res_verify = verify(sys_verify,params_verify,options,spec_verify);
tComp_verify = toc;
disp("Verification successful? " + res_verify);


% 2. falsification
sat = false;

[sys_falsify,params_falsify,spec_falsify,sat] = generateRandomVerify(...
    "StateDimension",n,...
    "InputDimension",m,...
    "OutputDimension",r,...
    "RealInterval",realInt,...
    "ImaginaryInterval",imagInt,...
    "NrSpecs",nrSpecs,...
    "SetRepInitialSet",setRepInitialSet,...
    "SetRepInputSet",setRepInputSet,...
    "SetRepSpec",setRepSpec,...
    "Satisfiable",sat, ...
    "RelError",relerror);

% write benchmark to JSON format
benchmark2json(sys_falsify,params_falsify,spec_falsify,'rand02.json');

% convert to zonotopes for verification algorithm
params_falsify.R0 = zonotope(params_falsify.R0);
params_falsify.U = zonotope(params_falsify.U);

% solve falsification task
tic;
res_falsify = verify(sys_falsify,params_falsify,options,spec_falsify);
tComp_falsify = toc;
disp("Falsification successful? " + ~res_falsify);

% output text
text{1} = ['Random,RAND01,',num2str(res_verify),',',num2str(tComp_verify)];
text{2} = ['Random,RAND02,',num2str(res_falsify),',',num2str(tComp_falsify)];

end

% ------------------------------ END OF CODE ------------------------------
