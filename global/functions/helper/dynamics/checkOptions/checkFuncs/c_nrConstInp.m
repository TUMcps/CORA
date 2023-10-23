function [res,msg] = c_nrConstInp(val,sys,params,options)
% c_nrConstInp - costum validation function to check whether
%    options.nrConstInp and params.u match
%
% Syntax:
%    [res,msg] = c_nrConstInp(val,sys,params,options)
%
% Inputs:
%    val - value for options.nrConstInp
%    sys - some contDynamics object
%    params - model parameters
%    options - algorithm parameters
%
% Outputs:
%    res - logical whether validation was successful
%    msg - error message if validation failed
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Authors:       Mark Wetzlinger
% Written:       10-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

% compute columns of u
u_cols = size(params.u,2);
if isa(sys,'linearSys') || isa(sys,'linearSysDT')
    if any(any(sys.D)) && size(params.u,2) > 1
        u_cols = size(params.u,2) - 1;
    end
end

isDT = isa(sys,'linearSysDT') || isa(sys,'nonlinearSysDT') || ...
                                                isa(sys,'neurNetContrSys');
if isDT
    reachSteps = round((params.tFinal - params.tStart)/sys.dt);
end

% check if size of params.u matches options.nrConstInp
if isscalar(val)
    if u_cols > 1
        % either equal to required length or integer multiple
        if mod(val,u_cols) ~= 0
            res = false;
            msg = ['does not comply with params.u: \n'...
                'If options.nrConstInp is a scalar value, and there is a\n'...
                'time-varying input vector params.u, it has to be\n'...
                'a multiple N*size(params.u,2) where N >= 1'];
            return;
            % same for u_cols == 1 -> combine!
        end
    elseif u_cols == 1
        % has to fit to sampling time: integer fraction / multiple
        if isDT && mod(reachSteps,val) ~= 0
            res = false;
            msg = ['does not comply with params.u: \n'...
                'For a discrete-time system, if options.nrConstInp is a scalar value\n'...
                'and there is no time-varying input params.u, it has to be a fraction steps/N,\n'...
                'where N >= 1 and steps = (tFinal-tStart)/dt'];
            return;
        end
        % for continuous-time and params.u only a vector, nrConstInp can be
        % an arbitrary integer value >= 1 (already checked)
    end
    % ...reqLength cannot be less than 1
else
    % options.nrConstInp is a vector: length has to match params.u (which
    % has already been ensured), values only have to be larger than 1
    % (which has been ensured in a previous check function)
    if u_cols ~= length(val)
        res = false;
        msg = ['does not comply with params.u: \n'...
            'The length of options.nrConstInp has to match the number of columns in params.u'];
        return;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
