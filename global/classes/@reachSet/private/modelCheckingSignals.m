function res = modelCheckingSignals(R,eq,varargin)
% modelCheckingSignals - check if a reachable set satisfies an STL formula
%                        using the three-valued signal approach
%
% Syntax:
%    res = modelCheckingSignals(R,eq)
%
% Inputs:
%    R - reachable set (class reachSet)
%    eq - logic formula (class stl)
%    varargin - Additional options as name-value pairs:
%               'returnBool' - return Boolean instead of Kleene result (default: true)
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
%    res = modelChecking(R,eq,'tvl')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Benedikt Seidl
% Written:       03-May-2023
% Last update:   08-February-2024 (FL, use stlInterval instead of interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check and extract additional options
NVpairs = varargin;
checkNameValuePairs(NVpairs, {'returnBool'});
[~, returnBool] = readNameValuePair(NVpairs, 'returnBool', 'islogical', true);

% prepare the formula
[phi,aps] = combineAtomicPropositions(desugar(disjunctiveNormalForm(eq)));

% compute the masks of this formula
msk = masks(phi,stlInterval(0,1),'any');

% get the duration of the reachable set
tFinal = query(R, 'tFinal');

% compute duration of resulting signal
dur = tFinal - maximumTime(phi);

% compute the signals for all atomic propositions
analyzer = reachSetAnalyzer(aps, tFinal, msk);
signals = analyzer.analyze(R);

% compute the validity of every path through the reachable set
res = kleene.True;

for i=1:length(signals)
    sig = evaluateSignal(phi,dur,signals{i});

    res = res & sig.at(0);
end

if returnBool
    res = res == kleene.True;
end

end

% ------------------------------ END OF CODE ------------------------------
