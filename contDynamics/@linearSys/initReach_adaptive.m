function [Rend,options,linOptions] = initReach_adaptive(sys,Rstart,options,linOptions)
% initReach_adaptive - computes the linearized reachable set
%  from an originally nonlinear system (linearization error = 0),
%  loops until time step found which satisfies abstraction error bound
%
% Syntax:
%    Rend = initReach_adaptive(sys,Rstart,linOptions)
%
% Inputs:
%    sys - linearized system
%    Rstart - reachable set of current time point
%    linOptions - options for linearized system
%    R - actual reachable set
%
% Outputs:
%    Rend - reachable set (time interval/point) of current time + time step
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%   -
%
% See also: linReach

% Authors:       Mark Wetzlinger
% Written:       25-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% exponential matrix and time interval error (incl. adaptive taylorTerms)
[sys,linOptions] = expmtie_adaptive(sys,linOptions);
% compute reachable set due to input
[sys,linOptions] = inputSolution(sys,linOptions);
%change the time step
sys.taylor.timeStep = linOptions.timeStep;
% save taylorTerms for analysis
options.tt_lin(options.i,1) = linOptions.taylorTerms;
options.etalinFro(options.i,1) = linOptions.etalinFro;

%compute reachable set of first time interval
eAt = expm(sys.A*linOptions.timeStep);
%save data to object structure
sys.taylor.eAt=eAt;

F=sys.taylor.F;
inputCorr=sys.taylor.inputCorr;
% if iscell(sys.taylor.Rtrans)
%     Rtrans=sys.taylor.Rtrans{1};
% else
    Rtrans=sys.taylor.Rtrans;
% end

%first time step homogeneous solution
Rhom_tp = eAt*Rstart + Rtrans;
Rhom = enclose(Rstart,Rhom_tp) + F*Rstart + inputCorr;

% preliminary solutions without RV
if isfield(options,'gredIdx')
    if length(options.gredIdx.Rhomti) == options.i
        % select the same generators for reduction as in the previous
        % iteration of the current time step to limit the influence of
        % the order reduction on the computation of the gain order
        Rend.ti = reduce(Rhom,'idx',options.gredIdx.Rhomti{options.i});
    else
        % reduce adaptively and store indices of generators that have
        % not been selected for reduction
        [Rend.ti,~,options.gredIdx.Rhomti{options.i}] = ...
            reduce(Rhom,'adaptive',options.redFactor);
    end
else
    % safety branch (call from algorithm without gredIdx)
    Rend.ti = reduce(Rhom,'adaptive',options.redFactor);
end

if isfield(options,'gredIdx')
    if length(options.gredIdx.Rhomtp) == options.i
        % select the same generators for reduction as in the previous
        % iteration of the current time step to limit the influence of
        % the order reduction on the computation of the gain order
        Rend.tp = reduce(Rhom_tp,'idx',options.gredIdx.Rhomtp{options.i});
    else
        % reduce adaptively and store indices of generators that have
        % not been selected for reduction
        [Rend.tp,~,options.gredIdx.Rhomtp{options.i}] = ...
            reduce(Rhom_tp,'adaptive',options.redFactor);
    end
else
    % safety branch (call from algorithm without gredIdx)
    Rend.tp = reduce(Rhom_tp,'adaptive',options.redFactor);
end

% reduce and add RV only if exists
if linOptions.isRV
    % read and reduce RV from struct
    if isfield(options,'gredIdx')
        if length(options.gredIdx.Rpar) == options.i
            % select the same generators for reduction as in the previous
            % iteration of the current time step to limit the influence of
            % the order reduction on the computation of the gain order
            RV = reduce(sys.taylor.RV,'idx',options.gredIdx.Rpar{options.i});
        else
            % reduce adaptively and store indices of generators that have
            % not been selected for reduction
            [RV,~,options.gredIdx.Rpar{options.i}] = ...
                reduce(sys.taylor.RV,'adaptive',options.redFactor);
        end
    else
        % safety branch (call from algorithm without gredIdx)
        RV = reduce(sys.taylor.RV,'adaptive',options.redFactor);
    end

    %total solution
    Rend.ti = Rend.ti + RV;
    Rend.tp = Rend.tp + RV;
end


end

% ------------------------------ END OF CODE ------------------------------
