function res = modelCheckingIncremental(R,phi,varargin)
% modelCheckingIncremental - incremental STL verification based on four-valued logic
%
% Syntax:
%    res = modelCheckingIncremental(R,phi)
%
% Inputs:
%    R - reachable set
%    phi - STL formula to verify
%    varargin - Additional options as name-value pairs:
%              'returnBool' - return Boolean instead of a four-valued result (default: true)
%              'propFreq' - number of observations to accumulate before propagating (default: 20)
%              'verbose' - print additional information (default: false)
%
% Outputs:
%    res - Boolean indicating whether the formula is satisfied
%         if res is false, the formula is not necessarily falsified,
%         the reachable set could also not be sufficient to verify the formula
%
% Example:
%    x = stl('x',1);
%    phi = finally(x(1) > 1,stlInterval(0,5));
%
%    sys = linearSys(0,1);
%
%    params.R0 = zonotope([0,0.5*eye(1)]);
%    params.U = zonotope([1,0.1*eye(1)]);
%    params.tFinal = 2;
%
%    options.timeStep = 0.01;
%    options.taylorTerms = 10;
%    options.zonotopeOrder = 20;
%
%    R = reach(sys,params,options);
%
%    res = modelChecking(R,phi,'incremental','propFreq',10,'verbose',true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl, incrementalMultiBranch, incrementalSingleBranch

% Authors:       Florian Lercher
% Written:       15-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check and extract additional options
NVpairs = varargin;
checkNameValuePairs(NVpairs,{'returnBool','propFreq','verbose'});
[NVpairs,returnBool] = readNameValuePair(NVpairs,'returnBool','islogical',true);
[NVpairs,propFreq] = readNameValuePair(NVpairs,'propFreq','isnumeric',20);
[~,verbose] = readNameValuePair(NVpairs,'verbose','islogical',false);

% create online analyzer
analyzer = onlineReachSetAnalyzer(desugar(phi),propFreq);

% check if the reachable set has multiple branches
if length(R) > 1
    tFinal = query(R,'tFinal');
    res = incrementalMultiBranch(R,analyzer,1,tFinal,verbose);
else
    res = incrementalSingleBranch(R,analyzer,verbose);
end

% convert to Boolean if necessary
if returnBool
    res = res == fourValued.True;
end

% ------------------------------ END OF CODE ------------------------------
