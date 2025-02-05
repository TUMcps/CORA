function res = modelChecking(R,eq,varargin)
% modelChecking - check if a reachable set satisfies an STL formula
%
% Syntax:
%    res = modelChecking(R,eq)
%    res = modelChecking(R,eq,alg)
%
% Inputs:
%    R - reachable set (class reachSet)
%    eq - logic formula (class stl)
%    alg - algorithm used, 'sampledTime' (Section 4.2 in [1]),
%          'rtl' (Theorem 1 in [1]), 'signals', or 'incremental'
%          (Corollary 1 in [2]). Note that some algorithms may accept
%          additional parameters.
%
% Outputs:
%    res - formula satisfied (true) or not (false)
%
% Example: 
%    x = stl('x',2);
%    eq = until(x(2) < -0.7,x(1) > 0.7,interval(0,2));
%    
%    sys = linearSys([0 -1; 1 0],[0;0]);
%
%    params.R0 = zonotope([0;-1]);
%    params.tFinal = 2;
%
%    options.timeStep = 0.5;
%    options.zonotopeOrder = 10;
%    options.taylorTerms = 10;
%
%    R = reach(sys,params,options);
%
%    res = modelChecking(R,eq)
%
% References: 
%    [1] H. Roehm et al. "STL Model Checking of Continuous and Hybrid
%        Systems", International Symposium on Automated Technology for 
%        Verification and Analysis, pp. 412-427, 2016.
%    [2] F. Lercher and M. Althoff, "Using Four-Valued Signal Temporal
%        Logic for Incremental Verification of Hybrid Systems", Computer
%        Aided Verification, 2024.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl, Florian Lercher
% Written:       09-November-2022 
% Last update:   15-February-2024 (FL, add incremental algorithm)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default algorithm
alg = setDefaultValues({'sampledTime'},varargin);

% check input arguments
inputArgsCheck({{R,'att','reachSet'}; ...
               {eq,'att','stl'};...
               {alg,'str',{'sampledTime','rtl','signals','incremental'}}});

% call the selected model checking algorithm
switch alg
    case 'sampledTime'
        res = priv_modelCheckingSampledTime(R,eq);
    case 'rtl'
        res = priv_modelCheckingRTL(R,eq);
    case 'signals'
        res = priv_modelCheckingSignals(R,eq,varargin{2:end});
    case 'incremental'
        res = priv_modelCheckingIncremental(R,eq,varargin{2:end});
end
    
% ------------------------------ END OF CODE ------------------------------
