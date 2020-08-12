function [Rout,Rout_tp,res,tVec] = reach_adaptive(obj,options)
% reach_adaptive - computes the reachable set for linear systems using the
%  propagation from the start as well as adaptive parametrization
%
% Syntax:
%    [Rout,Rout_tp,res,tVec] = reach_adaptive(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rout - output set of time intervals for the continuous dynamics 
%    Rout_tp - output set of points in time for the continuous dynamics
%    res - boolean (only if specification given)
%    tVec - vector containing time step sizes (e.g. for plotting)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       23-July-2019
% Last update:   08-Oct-2019
% Last revision: ---


%------------- BEGIN CODE --------------

% setup for calculation of Rout(_tp) --------------------------------------
[C,D,k] = initOutputEquation(obj,options);
% -------------------------------------------------------------------------

% initialize parameters ---------------------------------------------------
options = rA_initTimeStepTaylorTerms(obj,options);
options.zonotopeOrder = Inf;
formerTSTT = [options.timeStep; options.taylorTerms];
tVec = options.timeStep;
% -------------------------------------------------------------------------

% init step ---------------------------------------------------------------
[R_0, options] = rA_initReach(obj, options);
% compute output set
[Rout{1},Rout_tp{1}] = outputSet(C,D,k,R_0,options);
% safety property check
if isfield(options,'specification')
    if ~check(options.specification,Rout{1})
        % violation
        res = false;
        return
    end
end
% -------------------------------------------------------------------------

% access to values for loop -----------------------------------------------
Rhom_0 = options.Rhom;
Rhom_tp_0 = options.Rhom_tp;
eADelta = options.Q;
% -------------------------------------------------------------------------

% main loop: post ---------------------------------------------------------
cycle = 20;
i = 1; % init iteration counter
while options.tFinal - options.t > eps
    % update iteration counter
    i = i + 1;
    if isfield(options,'verbose') && options.verbose
        if mod(i,cycle) == 0
            fprintf("Step " + i + ": " + options.t + "\n");
        end
    end
    
    % find new time step / Taylor terms -----------------------------------
    % updated: options.timeStep|taylorTerms|t
    %          obj.taylor.powers|E|F
    options = rA_postTimeStepTaylorTerms(obj, options);
    tVec(i,1) = options.timeStep;
    % ---------------------------------------------------------------------
    
    
    % calculate new init procedure for current step -----------------------
    if any([options.timeStep; options.taylorTerms] ~= formerTSTT)
        [options, eADelta] = rA_post(obj, options);
        % access values for propagation
        Rhom_0 = options.Rhom;
        Rhom_tp_0 = options.Rhom_tp;
        formerTSTT = [options.timeStep; options.taylorTerms];
    end
    % ---------------------------------------------------------------------
    
    
    % propagate homogeneous part ------------------------------------------
    % H^R([t_i, t_i + Delta t_i]) = e^At_i * H^R([0, Delta t_i])
    Rhom    = options.Q * Rhom_0;
    % H^R([0, Delta t_i) = e^At_i * H^R(Delta t_i)
    Rhom_tp = options.Q * Rhom_tp_0;
    % ---------------------------------------------------------------------
    
    % propagate inhomogeneous part ----------------------------------------
    options = rA_propInhom(obj, options);
    % ---------------------------------------------------------------------
    
    % calculate full reachable set ----------------------------------------
    % R([t_i, t_i + Delta t_i]) = H([t_i, t_i + Delta t_i]) + P([0, t_i])
    Rcont.ti = Rhom + options.Rinhom;
    % R(t_i + Delta t_i) = H(t_i + Delta t_i) + P([0, t_i])
    Rcont.tp = Rhom_tp + options.Rinhom;
    % ---------------------------------------------------------------------
    
    
    % output set and safety property check --------------------------------
    [Rout{i,1},Rout_tp{i,1}] = outputSet(C,D,k,Rcont,options);
    % safety property check
    if isfield(options,'specification')
        if ~check(options.specification,Rout{i})
            % violation
            Rout = Rout(1:i);
            Rout_tp = Rout_tp(1:i);
            res = false;
            return
        end
    end
    % ---------------------------------------------------------------------
    
    
    % propagate matrix exponentials ---------------------------------------
    options.P = options.Q;
    options.Q = options.Q * eADelta;
    % ---------------------------------------------------------------------
    
end


% specification fulfilled at all times
res = true;

end


%------------- END OF CODE --------------