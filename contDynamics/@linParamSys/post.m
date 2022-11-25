function [Rnext,options] = post(obj,R,options)
% post - computes the reachable continuous set for one time step of a
% linear interval system
%
% Syntax:  
%    [Rnext] = post(obj,R,options)
%
% Inputs:
%    obj - linIntSys object
%    R - reachable set of the previous time step
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      15-May-2007 
% Last update:  07-January-2009
%               22-June-2009
%               29-June-2009
%               06-August-2010
%               08-August-2016
%               09-February-2016
%               19-May-2020 (MW, error handling for exploding sets)
% Last revision: ---

%------------- BEGIN CODE --------------

% rename for readability
M1 = obj.mappingMatrixSet.zono;
M2 = obj.mappingMatrixSet.int;

% for explosion detection (end of function)
if isa(R.ti,'zonotope')
    sizeRti = sum(abs(generators(R.ti)),2);
end

% use different algorithm if the input changes in every time step
if isfield(options,'uTransVec')
   
    % update input sets
    inputSolution(obj,options);
    
    % compute time point solution
    Rinit = R.tp;
    Rhom_tp = obj.mappingMatrixSet.zono*Rinit + obj.mappingMatrixSet.int*Rinit;
    
    % compute time interval solution
    if isa(Rinit,'zonoBundle')
        Rhom = enclose(Rinit,Rhom_tp+obj.Rtrans) ...
                + obj.F*Rinit.Z{1} + obj.inputCorr + (-1*obj.Rtrans);
    else
        Rhom = enclose(Rinit,Rhom_tp+obj.Rtrans) ...
                + obj.F*Rinit + obj.inputCorr + (-1*obj.Rtrans);
    end

    % total solution
    Rnext.ti = reduce(Rhom + obj.Rinput,...
        options.reductionTechnique,options.zonotopeOrder);
    Rnext.tp = reduce(Rhom_tp + obj.Rinput,...
        options.reductionTechnique,options.zonotopeOrder);
    
else

    % next set (no input)
    R.ti = M1*R.ti + M2*R.ti; 

    % bloating due to input
    R.ti = R.ti + obj.Rinput;

    % reduce zonotope
    Rred=reduce(R.ti,options.reductionTechnique,options.zonotopeOrder);

    % write results to reachable set struct Rnext
    Rnext.ti=Rred;

    if options.compTimePoint
        % next set (no input)
        R.tp = M1*R.tp + M2*R.tp;

        % bloating due to input
        R.tp = R.tp + obj.Rinput;

        % reduce zonotope
        Rred_tp=reduce(R.tp,options.reductionTechnique,options.zonotopeOrder);

        % write results to reachable set struct Rnext
        Rnext.tp=Rred_tp;
    else
        Rnext.tp = [];
    end
end

% check for explosion
if isa(R.ti,'zonotope')
    if max( sum(abs(generators(Rnext.ti)),2) ./ sizeRti ) > 100 % arbitary value
        throw(printExplosionError());
    end
end


%------------- END OF CODE --------------