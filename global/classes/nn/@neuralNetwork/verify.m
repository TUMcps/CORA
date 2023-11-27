function [res, x] = verify(nn, X0, spec, varargin)
% verify - automated verifciation for specification on neural networks.
%    Note: If the input set is not an interval, a potentially found
%    falsifying example might not be in the actual input set.
%    However, if the verification is successful, also the actual input set
%    is verified.
%
% Syntax:
%    [res, x] = verify(nn, X0, spec)
%    [res, x] = verify(nn, X0, spec, varargin)
%
% Inputs:
%    nn - object of class neuralNetwork
%    X0 - initial set of class interval
%    spec - specifications of class specification
%    varargin - (optional) name-value pairs
%       - 'Splits': maximum number of recursive splits of X0 
%       - 'RefinementSteps': number of refinements per step
%       - 'Verbose': true/false, verbose log
%       - 'Plot': true/false, plot verification progress
%                   (will increase the verification time)
%       - 'PlotDimsIn': dimensions to plot in input space
%       - 'PlotDimsOut': dimensions to plot in output space
%
% Outputs:
%    res - result: true if specification is satisfied, false if not, empty if unknown
%    x - counterexample in terms of an initial point violating the specs
%
% References:
%    [1] Ladner, T., et al. (2023). Automatic abstraction refinement in
%        neural network verification using sensitivity analysis. HSCC '23:
%        Proceedings of the 26th International Conference on
%        Hybrid Systems: Computation and Control.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       23-November-2021
% Last update:   30-November-2022 (TL, removed neuralNetworkOld, adaptive)
%                25-July-2023 (TL, input parsing, improvements)
%                23-November-2023 (TL, verbose, bug fix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
if nargin < 3
    throw(CORAerror("CORA:notEnoughInputArgs", 3))
end
[nn,X0,spec,splits,numRefineSteps,verbose,doPlot,plotDimsIn,plotDimsOut] = ...
    aux_parseInput(nn,X0,spec,varargin{:});

% 2. preprocess
[nn,X0,spec] = aux_preprocess(nn,X0,spec);

% 3. initial plotting
if doPlot
    [xLim, yLim, hanSpec] = aux_plotInit(nn,X0,spec,plotDimsIn,plotDimsOut);
end

% 4. start evaluation -----------------------------------------------------

% init unknown sets
unknownSets = {X0};

% main evaluation loop
if verbose
    disp("Verifying neural network ...")
end
while ~isempty(unknownSets) && splits >= 0
    if verbose
        fprintf("Nr. of remaining splits: %2d - Nr. of unknown sets: %3d\n", splits, length(unknownSets))
    end
    
    splits = splits-1;

    % iterate over all unknown sets
    newUnknownSets = cell(1, 2*length(unknownSets));
    unknownCnt = 0;
    for i = 1:length(unknownSets)
        X = unknownSets{i};

        % a) try to falsify set -------------------------------------------
        [res, x] = aux_falsify(nn, X, spec);
        if ~res 
            % found falsifying point from input set
            if doPlot
                markerSize = 50;
                subplot(1, 2, 1);
                scatter(x(plotDimsIn(1)), x(plotDimsIn(2)), markerSize, 'oc', 'filled', 'DisplayName', 'Counterexample')
                subplot(1, 2, 2);
                y = nn.evaluate(x);
                scatter(y(plotDimsOut(1)), y(plotDimsOut(2)), markerSize, 'oc', 'filled', 'DisplayName', 'Counterexample')
                drawnow
            end

            if verbose
                disp("Result: FALSIFIED")
            end
            return;
        end

        % b) try to verify set --------------------------------------------
        [res, Y] = aux_verify(nn, X, spec, numRefineSteps);

        if res
            % set already verified, no further processing required

            if doPlot
                % plot verified sets
                subplot(1, 2, 1);
                plot(X, plotDimsIn, 'Color', colorblind('y'), 'HandleVisibility', 'off');
                subplot(1, 2, 2);
                plot(Y, plotDimsOut, 'Color', colorblind('y'), 'HandleVisibility', 'off');

                % update spec plot if limits changed
                if ~all([xLim yLim] == [xlim() ylim()])
                    % replot specifications underneath all other plots
                    hanSpecOld = hanSpec;
                    name = hanSpecOld.DisplayName;
                    hanSpec = plot(spec, plotDimsOut, 'DisplayName',name);
                    uistack(hanSpec,'bottom')
                    % delete old specification (after plotting spec to
                    % avoid axis rescaling)
                    delete(hanSpecOld)

                    % update limits
                    xLim = xlim;
                    yLim = ylim;
                end

                % drawnow
                drawnow
            end

            continue;
        end

        % set not verified

        % if doPlot
        %     % plot unverified set
        %     subplot(1, 2, 1);
        %     plot(X, plotDimsIn, 'Color', colorblind('b'), 'HandleVisibility', 'off');
        %     subplot(1, 2, 2);
        %     plot(Y, plotDimsOut, 'Color', colorblind('b'), 'HandleVisibility', 'off');
        %     drawnow
        % end

        % c) split X to get a tighter over-approximation ------------------
        [X1, X2] = aux_split(nn, X, spec, Y);
        newUnknownSets(unknownCnt + [1,2]) = {X1, X2};
        unknownCnt = unknownCnt + 2;

        if doPlot
            subplot(1, 2, 1);
            plot(X1, plotDimsIn, 'Color', colorblind('b'), 'HandleVisibility', 'off');
            plot(X2, plotDimsIn, 'Color', colorblind('b'), 'HandleVisibility', 'off');
            drawnow
        end

    end

    % update unknown sets with unknown splitted sets
    unknownSets = newUnknownSets(1:unknownCnt);
end

if isempty(unknownSets)
    % input set X0 verified
    if verbose
        disp("Result: VERIFIED")
    end
    res = true;
else
    % unknown result
    if verbose
        disp("Result: UNKNOWN")
    end
    res = [];
end

end


% Auxiliary functions -----------------------------------------------------

function [nn,X0,spec,splits,numRefineSteps,verbose,doPlot,plotDimsIn,plotDimsOut] = ...
    aux_parseInput(nn,X0,spec,varargin)

    % parse name-value pairs
    [varargin,splits] = readNameValuePair(varargin,'Splits','isscalar',10);
    [varargin,numRefineSteps] = readNameValuePair(varargin,'RefinementSteps',{'isnumeric','isscalar'},1);
    [varargin,verbose] = readNameValuePair(varargin,'Verbose',{'islogical','isscalar'},false);
    [varargin,doPlot] = readNameValuePair(varargin,'Plot',{'islogical','isscalar'},false);
    [varargin,plotDimsIn] = readNameValuePair(varargin,'PlotDimsIn','isnumeric',[1 2]);
    [~,plotDimsOut] = readNameValuePair(varargin,'PlotDimsOut','isnumeric',[1 2]);

    if CHECKS_ENABLED
        inputArgsCheck({; ...
            {nn, 'att', {'neuralNetwork'}}; ...
            {X0, 'att', {'contSet'}}; ...
            {spec, 'att', {'specification'}}; ...
            });

        % check input set
        if dim(X0) ~= nn.neurons_in
            throw(CORAerror("CORA:wrongValue",'second','Number of input neurons and dimension of input set must match.'))
        end
        
        % check specification
        for i=1:length(spec)
            if ~ismember(spec(i).type, {'safeSet','unsafeSet'})
                throw(CORAerror("CORA:wrongValue",'third','Only safe and unsafe specifications are supported.'))
            end

            if dim(spec(i).set) ~= nn.neurons_out
                throw(CORAerror("CORA:wrongValue",'third','Number of output neurons and dimension of specification sets must match.'))
            end

            if ~representsa_(spec(i).time,'emptySet',eps) 
                throw(CORAerror("CORA:wrongValue",'third','Time of specifications will be ignored.'))
            end
        end

        % parse name-value pairs
        if length(plotDimsIn) ~= 2
            throw(CORAerror("CORA:wrongValue",'name-value pair ''PlotDimsIn''','Must specify 2 dimensions.'))
        end
        if length(plotDimsOut) ~= 2
            throw(CORAerror("CORA:wrongValue",'name-value pair ''PlotDimsOut''','Must specify 2 dimensions.'))
        end
    end
end

function [nn,X0,spec] = aux_preprocess(nn,X0,spec)

    % transform input to an interval
    if ~representsa_(X0,'interval',eps)
        warning("CORA Neural Network Verification: "+ ...
            "Input Set is not an interval. "+ ...
            "Thus, a potentially found falsifying example "+ ...
            "might not be in the actual input set. "+ ...
            "However, if the verification is successful, "+ ...
            "also the actual input set is verified.")
    end    
    X0 = interval(X0);

    % classification task can be significantly accelerated
    % by applying the argmax trick [1, Prop. B.2]
    % -> check if it is a classification task
    if length(spec) == 1 && isa(spec.set, 'polytope')

        % read polytope set
        P = spec.set;
        
        % use P.A as final linear layer of the neural network
        nn = neuralNetwork([ ...
            nn.layers(:); {nnLinearLayer(P.A, 0, 'Argmax Trick')}
        ]);

        spec = specification(interval(-inf(length(P.b),1),P.b),spec.type);
    end
    
    % otherwise keep as is
end

function [xLim, yLim, hanSpec] = aux_plotInit(nn,X0,spec,plotDimsIn,plotDimsOut)

    figure;
    subplot(1, 2, 1);
    hold on;
    title("Input")
    plot(X0, plotDimsIn, '-', 'Color', colorblind('b'), 'HandleVisibility', 'off');

    xs = [X0.randPoint(1000), X0.randPoint(100, 'extreme')];
    ys = nn.evaluate(xs);

    plot(xs(plotDimsIn(1), 1), xs(plotDimsIn(2), 1), '-', 'Color', colorblind('b'), 'DisplayName', 'Unknown')
    plot(xs(plotDimsIn(1), 1), xs(plotDimsIn(2), 1), '-', 'Color', colorblind('y'), 'DisplayName', 'Verified')
    scatter(xs(plotDimsIn(1), :), xs(plotDimsIn(2), :), '.k', 'DisplayName', 'Samples');
    legend();

    subplot(1, 2, 2);
    hold on; legend();
    title("Output")
    plot(ys(plotDimsOut(1), 1), ys(plotDimsOut(2), 1), '-', 'Color', colorblind('y'), 'DisplayName', 'Verified')
    scatter(ys(plotDimsOut(1), :), ys(plotDimsOut(2), :), '.k', 'DisplayName', 'Samples');

    % retrieve limits (before plotting specs!)
    xLim = xlim;
    yLim = xlim;

    % plot spec at last for correct axis limits; 
    % move underneath other plots though
    if strcmp(spec(1).type, 'unsafeSet')
        specName = 'Unsafe Specification';
    else
        specName = 'Safe Specification';
    end
    hanSpec = plot(spec, plotDimsOut, 'DisplayName', specName);
    uistack(hanSpec,'bottom')

    % drawnow
    drawnow
end

function [res, x] = aux_falsify(nn, X, spec)
% try to falsify the specifications

% sample pounts % TODO critical point
xs = [X.randPoint(1000), X.randPoint(100, 'extreme')];
ys = nn.evaluate(xs);

% test points in output space
[res,~,indObj] = check(spec, ys);
if res
    % res = true;
    x = [];

else
    % res = false; get counter example
    x = xs(:, indObj);

end

end

function [res, Yout] = aux_verify(nn, X, spec, numRefineSteps)
% try to verify the specifications

% TODO: more advanced methods? Niklas had in neuralNetworkOld a more
% advanced version depending on the specific set representation used
% for X and spec. Might be useful to look into. (see remove commit).

% evaluation parameters
evParams = struct;
% evParams.reuse_bounds = true;
evParams.num_generators = 10000;
evParams.reuse_bounds = true;

nn.reset();
% evaluate neural network
for i = 1:numRefineSteps
    X = polyZonotope(X);
    Y = nn.evaluate(X, evParams);

    if i == 1
        % return linearized over-approximation if unable to verify
        Yout = Y;
    end

    if isa(spec(1).set, 'interval')
        Y = interval(Y);
    end

    % check specification
    res = check(spec, Y);
    if res % verified successfully
        Yout = Y;
        return
    end

    % verify
    nn.refine(2, 'layer', 'both', X.randPoint(1));
end

% unable to verify
res = false;
end

function [X1, X2] = aux_split(nn, X, spec, Y)
% splits the given set

% TODO: more advanced methods? Niklas had in neuralNetworkOld a more
% advanced version depending on the specific set representation used
% for X and spec. Might be useful to look into. (see remove commit).

type = 'sensitivity';
if strcmp(type, 'longest')
    % split at longest axis
    [~, n] = max(X.sup-X.inf);

elseif strcmp(type, 'sensitivity')
    % split at axis with largest sensitivity
    x = X.randPoint(1);
    S = nn.calcSensitivity(x);
    S = vecnorm(S, 2, 1)' + eps; % eps to avoid 0 sensitivity
    [~, n] = max(S .* (X.sup - X.inf)); % multiply with axis length
else
    throw(CORAerror('CORA:wrongValue', sprintf("Unknown splitting type: %s.", type)))
end

Xs = split(X, n);
X1 = Xs{1};
X2 = Xs{2};
end

% ------------------------------ END OF CODE ------------------------------
