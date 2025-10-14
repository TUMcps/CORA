function res = monitorSTL(traj,eq)
% monitorSTL - check if a simulation result satisfies an STL formula
%
% Syntax:
%    res = monitorSTL(traj,eq)
%
% Inputs:
%    traj - simulation result (class trajectory)
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
%    traj = simulateRandom(sys, params, options);
%
%    res = monitorSTL(traj,eq)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Benedikt Seidl, Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{traj,'att',{'trajectory'}};...
                {eq,'att',{'stl'}}});

[phi,aps] = combineAtomicPropositions(desugar(eq));

res = true;


for i = 1:length(traj)
    analyzer = trajectoryAnalyzer(aps,traj(i));
    for s = 1:size(traj(i).x,3)
        sigs = analyzer.analyze(s);
        dur = traj(i).t(:,end);
    
        sig = evaluateSignal(phi,dur-maximumTime(phi),sigs,'logical');
    
        res = res && sig.at(0);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
