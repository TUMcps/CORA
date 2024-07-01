function [sys,params,spec,sat] = generateRandomVerify(varargin)
% generateRandomVerify - generates a random verification benchmark
%
% Syntax:
%    [sys,params,spec,sat] = generateRandomVerify(...
%                               "StateDimension",n,
%                               "InputDimension",m,
%                               "OutputDimension",r,
%                               "RealInterval",realInt,
%                               "ImaginaryInterval",imagInt,
%                               "ConvergenceFactor",convFactor,
%                               "NrSpecs",nrSpecs,
%                               "SetRepSpec",setRepSpec,
%                               "Satisfiable",sat,
%                               "MaxError",err)
%                       
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'StateDimension',n> - state dimension
%       <'InputDimension',nrInputs> - input dimension
%       <'OutputDimension',nrOutputs> - output dimension
%       <'RealInterval',realInt> - interval for real part of eigenvalues
%       <'ImaginaryInterval',imagInt> - interval for imaginary part of
%                                           eigenvalues
%       <'ConvergenceFactor',convFactor> - convergence factor for time
%                                           horizon
%       <'NrSpecs',nrSpecs> - number of specifications
%       <'SetRepSpec',setRepSpec> - convex set representation class
%       <'Satisfiable',sat> - true/false whether specifications are
%                               satisfiable
%       <'MaxError',maxerror> - difficulty of benchmark, i.e.,
%           tightness of reachable sets required to obtain correct result
%       <'RelError',relerror> - difficult of benchmark via error relative
%           to the size of the initial set
%
% Outputs:
%    sys - linearSys object
%    params - model parameters
%       .R0 - initial set
%       .U - input set
%       .tFinal - time horizon
%    spec - array of specifications
%    sat - true/false whether specifications are satisfiable
%
% Example: 
%    [sys,params,spec,sat] = generateRandomVerify("StateDimension",5);

% Authors:       Mark Wetzlinger
% Written:       06-May-2024
% Last update:   07-June-2024 (MW, use maximum/relative error for difficulty)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. read name-value pairs
% name-value pairs -> number of input arguments is always a multiple of 2
if nargin == 0
    throw(CORAerror('CORA:notEnoughInputArgs',2));
elseif mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
end
[n,nrInputs,nrOutputs,realInt,imagInt,convFactor,nrSpecs,...
    setRepInitialSet,setRepInputSet,setRepSpec,sat,maxError,relError] = ...
    aux_readNameValuePairs(varargin);

% 2. generate random linear system (use default values from there)
sys = linearSys.generateRandom(...
    'StateDimension',n,...
    'InputDimension',nrInputs,...
    'OutputDimension',nrOutputs,...
    'RealInterval',realInt,...
    'ImaginaryInterval',imagInt);

% 3. set parameters
[params,options] = aux_setParamsOptions(sys.A,nrInputs,...
    setRepInitialSet,setRepInputSet,convFactor,maxError,relError);

% 4. compute reachable set
R = aux_reach(sys,params,options);

% 5. place specifications
timeInt = interval(0,params.tFinal);
spec = aux_placeSpecs(R,sat,options.error,nrSpecs,setRepSpec,timeInt);

end


% Auxiliary functions -----------------------------------------------------

function [n,nrInputs,nrOutputs,realInt,imagInt,convFactor,nrSpecs,...
    setRepInitialSet,setRepInputSet,setRepSpec,...
    sat,maxerror,relerror] = aux_readNameValuePairs(NVpairs)

% check list of name-value pairs
checkNameValuePairs(NVpairs,{'StateDimension','InputDimension',...
    'OutputDimension','RealInterval','ImaginaryInterval',...
    'ConvergenceFactor','NrSpecs','SetRepInitialSet','SetRepInputSet',...
    'SetRepSpec','Satisfiable','MaxError','RelError'});

% state dimension given?
[NVpairs,n] = readNameValuePair(NVpairs,'StateDimension');

% input dimension given?
[NVpairs,nrInputs] = readNameValuePair(NVpairs,'InputDimension');

% output dimension given?
[NVpairs,nrOutputs] = readNameValuePair(NVpairs,'OutputDimension');

% interval for real part of eigenvalues given?
[NVpairs,realInt] = readNameValuePair(NVpairs,'RealInterval',...
    @(x) representsa_(x,'interval',eps),interval.empty(1));

% interval for imaginary part of eigenvalues given?
[NVpairs,imagInt] = readNameValuePair(NVpairs,'ImaginaryInterval',...
    @(x) representsa_(x,'interval',eps),interval.empty(1));

% convergence factor for time horizon given?
[NVpairs,convFactor] = readNameValuePair(NVpairs,'ConvergenceFactor',...
    @isnumeric,0.001);

% number of specifications given?
[NVpairs,nrSpecs] = readNameValuePair(NVpairs,'NrSpecs',...
    @isnumeric,1);

% set representation of initial set
admissibleContSetR0 = {'interval','zonotope'};
[NVpairs,setRepInitialSet] = readNameValuePair(NVpairs,'SetRepInitialSet',...
    @(x) any(strcmp(x,admissibleContSetR0)),'zonotope');

admissibleContSetU = {'interval','zonotope'};
[NVpairs,setRepInputSet] = readNameValuePair(NVpairs,'SetRepInputSet',...
    @(x) any(strcmp(x,admissibleContSetU)),'zonotope');

% set representation of specifications given?
admissibleContSetSpec = {'interval','zonotope','ellipsoid','halfspace'};
[NVpairs,setRepSpec] = readNameValuePair(NVpairs,'SetRepSpec',...
    @(x) any(strcmp(x,admissibleContSetSpec)),'interval');

% satisfiability prescribed?
[NVpairs,sat] = readNameValuePair(NVpairs,'Satisfiable',...
    @islogical,true);

% maximum Hausdorff distance prescribed?
[NVpairs,maxerror] = readNameValuePair(NVpairs,'MaxError',...
    @isnumeric,1);

% maximum relative Hausdorff distance prescribed?
[NVpairs,relerror] = readNameValuePair(NVpairs,'RelError',...
    @isnumeric,1);

end

function [params,options] = aux_setParamsOptions(A,nrInputs,...
    setRepInitialSet,setRepInputSet,convFactor,maxError,relError)

% read out dimension
n = length(A);

% sample initial set
switch setRepInitialSet
    case 'zonotope'
        params.R0 = zonotope.generateRandom("Dimension",n,"Center",5*ones(n,1));
    case 'interval'
        params.R0 = interval.generateRandom("Dimension",n,"Center",5*ones(n,1));
end
params.R0 = enlarge(params.R0,1./8*rad(interval(params.R0)));

% sample input set
switch setRepInputSet
    case 'zonotope'
        params.U = 0.01 * zonotope.generateRandom("Dimension",nrInputs);
    case 'interval'
        params.U = 0.01 * interval.generateRandom("Dimension",nrInputs);
end
params.U = enlarge(params.U,1/100*max(rad(interval(params.R0))));

% compute time horizon
params.tFinal = aux_computeTimeHorizon(A,convFactor);

% set Hausdorff error
if ~isempty(maxError)
    options.error = maxError;
else
    options.error = vecnorm(sum(2*rad(interval(params.R0)))) * relError;
end
% fix adaptive algorithm
options.linAlg = 'adaptive';

end

function tFinal = aux_computeTimeHorizon(A,convFactor)
% compute the time horizon for the randomly generated system using the
% largest eigenvalue of the system matrix and a convergence factor

% compute maximum real-part of eigenvalues of system matrix
ev = eigs(A);
ev_max = max(real(ev));

% time horizon until convergence
tFinal = log(convFactor) / ev_max;

end

function R = aux_reach(sys,params,options)

% convert to zonotope
params.R0 = zonotope(params.R0);
params.U = zonotope(params.U);

% compute outer approximation
R = reach(sys,params,options);

end

function spec = aux_placeSpecs(R,sat,maxError,nrSpecs,setRepSpec,timeInt)
% place the specifications

% read out time-interval solutions
Rlist = R.timeInterval.set;

% dimension
nrOutputs = dim(Rlist{1});
% list of directions (to avoid using similar directions twice)
ell_all = zeros(nrOutputs,nrSpecs);
% maximum tries per direction
max_iterations = nrSpecs + 100;
iteration = 1;
spec = [];

while length(spec) < nrSpecs && iteration <= max_iterations
    % force current iteration if not enough left to play around
    forceThisIteration = iteration >= (max_iterations - (nrSpecs - length(spec)));

    % choose random direction (unit length)
    ell = randn(nrOutputs,1);
    ell = ell ./ vecnorm(ell);

    ell_all(:,iteration) = ell;
    if iteration > 1 && any(ell' * ell_all(:,1:length(spec)) >= 0.9) && ...
            ~forceThisIteration
        % the closer the dot product is to 1, the more similar the
        % current direction is to a previously chosen one
        iteration = iteration + 1;
        continue
    end

    % compute support function and associated support vector of the
    % reachable set in the chosen direction        
    [~,sV_R,hitsInitSet] = aux_supportFunc(Rlist,ell);
    if hitsInitSet && ~forceThisIteration
        iteration = iteration + 1;
        continue
    end

    % random set centered at the origin
    switch setRepSpec
        case 'interval'
            S = interval.generateRandom("Dimension",nrOutputs,...
                "Center",zeros(nrOutputs,1));
        case 'zonotope'
            S = zonotope.generateRandom("Dimension",nrOutputs,...
                "Center",zeros(nrOutputs,1));
        case 'ellipsoid'
            S = ellipsoid.generateRandom("Dimension",nrOutputs,...
                "Center",zeros(nrOutputs,1));
        case 'polytope'
            S = polytope.generateRandom("Dimension",nrOutputs);
            S = S - center(S);
        case 'halfspace'
            % no randomness here...
            S = halfspace(-ell,0);
            sV_S = 0;
        otherwise
            throw(CORAerror('CORA:notSupported'));
    end

    % compute support vector in opposite direction
    if ~strcmp(setRepSpec,'halfspace')
        [~,sV_S] = supportFunc_(S,-ell,'upper');
    end

    % move the center of the specification so that is satisfiable
    S = S + (sV_R - sV_S);

    % for unsatisfiable specifications, we enlarge the specification by a
    % ball with the radius of the Hausdorff distance
    if ~sat
        % note: other than for halfspaces, this does not fulfill the
        % Hausdorff error anymore (i.e., the specifications are easier to
        % falsify than what the given error implies), as the set
        % representation are not closed under Minkowski sum with a ball
        switch setRepSpec
            case 'halfspace'
                S = halfspace(S.c, S.d + vecnorm(S.c*maxError));
            case 'zonotope'
                S = S + zonotope(zeros(nrOutputs,1),maxError*eye(nrOutputs));
            case 'interval'
                S = S + interval(-maxError*ones(nrOutputs,1),...
                    maxError*ones(nrOutputs,1));
            case 'ellipsoid'
                S = S + ellipsoid(maxError^2*eye(nrOutputs),zeros(nrOutputs,1));
            case 'polytope'
                S = polytope(S.A, S.b + vecnorm(S.A*maxError));
        end
    end

    % specifications should not intersect one another (unless for
    % halfspaces, where one cannot really avoid that)
    if ~strcmp(setRepSpec,'halfspace')
        intersectionFound = false;
        for i=1:length(spec)
            if isIntersecting_(S,spec(i).set,'approx')
                intersectionFound = true; break;
            end
        end
        % try again unless we would otherwise reach the end of the while
        % loop before generating the requested amount of specifications
        if intersectionFound && ~forceThisIteration
            continue
        end
    end

    % append to list of specifications
    if iteration == 1
        spec = specification(S,'unsafeSet',timeInt);
    else
        spec = [spec; specification(S,'unsafeSet',timeInt)];
    end

    % increment iteration counter
    iteration = iteration + 1;
end

end

function [sF,sV,hitsInitSet] = aux_supportFunc(Rlist,ell)
% compute the support function value and support vector of the reachable
% set (given as the union of individual sets) in a given direction;
% depending on the satisfiability of the subsequently generated
% specifications, we have to use the outer or inner approximation,
% accordingly

% read out number of sets, init support function value
nrSets = length(Rlist);
sF = -Inf;
idx = 0;

% loop over all reachable sets
for k=1:nrSets
    % compute support function of k-th set
    [sF_i,sV_i] = supportFunc_(Rlist{k},ell,'upper');

    % check value, override if larger
    if sF_i > sF
        idx = k;
        sF = sF_i;
        sV = sV_i;
    end
end

% is the maximum extent of the reachable set in the chosen direction given
% by the initial set?
hitsInitSet = idx == 1;

end

% ------------------------------ END OF CODE ------------------------------
