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
%    alg - algorithm used, 'sampledTime' (Section 4.2 in [1]) or 
%           'rtl' (Theorem 1 in [1]) or 'signals'
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
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default algorithm
alg = setDefaultValues({'sampledTime'},varargin);

% check input arguments
inputArgsCheck({{R,'att','reachSet'}; ...
               {eq,'att','stl'};...
               {alg,'str',{'sampledTime','rtl','signals'}}});

% call the selected model checking algorithm
if strcmp(alg,'sampledTime')
    res = modelCheckingSampledTime(R,eq);
elseif strcmp(alg,'rtl')
    res = modelCheckingRTL(R,eq);
elseif strcmp(alg,'signals')
    res = modelCheckingSignals(R,eq);
end
    
% ------------------------------ END OF CODE ------------------------------
