function res = monitorSTL(simRes,eq)
% monitorSTL - check if a simulation result satisfies an STL formula
%
% Syntax:
%    res = monitorSTL(simRes,eq)
%
% Inputs:
%    simRes - simulation result (class simResult)
%    eq - logic formula (class stl)
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
%    options.points = 5;
%
%    simRes = simulateRandom(sys, params, options);
%
%    res = monitorSTL(simRes,eq)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Benedikt Seidl
% Written:       14-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{simRes,'att',{'simResult'},{''}};...
                {eq,'att',{'stl'},{''}}});

[phi,aps] = combineAtomicPropositions(desugar(eq));

res = true;

analyzer = simResultAnalyzer(aps,simRes);

for i = 1:length(simRes.t)
    sigs = analyzer.analyze(i);
    dur = simRes.t{i}(end);

    sig = evaluateSignal(phi,dur-maximumTime(phi),sigs,'logical');

    res = res && sig.at(0);
end

end

% ------------------------------ END OF CODE ------------------------------
