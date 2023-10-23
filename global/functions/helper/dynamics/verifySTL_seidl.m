function res = verifySTL_seidl(sys,params,options,eq)
% verifySTL_seidl - fully automatically verify the given system
%
% Syntax:
%    res = verifySTL_seidl(sys,params,options,eq)
%
% Inputs:
%    sys - continuous or hybrid system
%    params - system parameters
%    options - reachability options
%    eq - STL formula
%
% Outputs:
%    res - formula satisfied (true) or not (false)
%
% Example:
%    x = stl('x',2);
%    eq = until(x(2) > 0.4,x(1) < 0,interval(0,2));
%
%    sys = linearSys([0 -1; 1 0],[0; 0]);
%
%    params.R0 = zonotope(interval([-1;-1],[1;1]));
%    options.taylorTerms = 10;
%    options.zonotopeOrder = 10;
%    options.verifyAlg = 'stl:seidl';
%
%    verify(sys,params,options,eq)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics/verify, hybridDynamics/verify, stl

% Authors:       Benedikt Seidl
% Written:       15-May-2023
% Last update:   ---
% Last revision: 19-October-2023 (TL, rename to verifySTL)

% ------------------------------ BEGIN CODE -------------------------------

% prepare formula ---------------------------------------------------------

% collect atomic propositions
[phi,aps] = combineAtomicPropositions(desugar(disjunctiveNormalForm(eq)));

% calculate duration and initial step size
dur = maximumTime(phi);

if ~isfield(options, 'timeStep')
    options.timeStep = dur / 2;
end

% compute the masks of this formula
msk = masks(phi,interval(0,options.timeStep),'any');

% prepare refinement loop -------------------------------------------------

params.tFinal = dur + options.timeStep;

% max 5 rounds
for i = 1:5
    [res,signals] = aux_reachAndVerify();

    % return if system is definitely valid
    if isequal(res, kleene.True)
        res = true;
        return
    end

    % return if system is definitely broken
    if isequal(res, kleene.False)
        res = false;
        return
    end

    % repeat reachability analysis with more accuracy when unknown
    options.timeStep = options.timeStep / 2;

    aux_updateMasks(signals);
end

res = false;


% Auxiliary functions -----------------------------------------------------

function [res,signals] = aux_reachAndVerify()

    % compute the reachable set
    R = reach(sys,params,options);

    % compute the signals for all atomic propositions
    analyzer = reachSetAnalyzer(aps,params.tFinal,msk);
    signals = analyzer.analyze(R);

    % compute the validity of every path through the reachable set
    res = kleene.True;

    for j=1:length(signals)
        sig = evaluateSignal(phi,params.tFinal - dur,signals{j});

        res = res & sig.at(0);
    end

end

function aux_updateMasks(signals)

    k = keys(msk);

    for j = 1:length(signals)
        sigs = signals{j};

        for l = 1:length(k)
            msk(k{l}) = combine(msk(k{l}),sigs(k{l}), @aux_comb);
        end
    end

end

function out = aux_comb(m,s)

    % only check previously unknown areas
    out = m & isequal(s,kleene.Unknown);

end

end

% ------------------------------ END OF CODE ------------------------------
