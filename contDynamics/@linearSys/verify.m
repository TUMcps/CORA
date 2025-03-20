function [res,varargout] = verify(linsys,params,options,spec)
% verify - verifies a hybrid system against the given specifications
%
% Syntax:
%    R = verify(linsys,params,options,spec)
%
% Inputs:
%    linsys - contDynamics object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%       .verifyAlg: 'reachavoid:zonotope','reachavoid:supportFunc',
%                   'stl:kochdumper','stl:seidl'
%    spec - object of class specification (reach-avoid) or stl
%
% Outputs:
%    res - true/false whether specifications are satisfied
%    varargout - depending on selected verification algorithm
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics/verify

% Authors:       Tobias Ladner
% Written:       19-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. check number of inputs
narginchk(4,4);

% 2. validate inputs
inputArgsCheck({ ...
    {linsys, 'att', 'linearSys'}; ...
    {params, 'att', 'struct'}; ...
    {options, 'att', 'struct'}; ...
    {spec, 'att', {'specification','stl'}}; ...
});

% 3. select verification algorithm
[verifyAlg,options] = aux_selectVerifyAlgorithm(linsys,params,options,spec);

% 5. call verify function depending on given specification
switch verifyAlg

    % reach-avoid
    case 'reachavoid:zonotope'
        varargout = cell(1,1);
        [res,varargout{:}] = priv_verifyRA_zonotope(linsys,params,options,spec);
    case 'reachavoid:supportFunc'
        varargout = cell(1,2);
        [res,varargout{:}] = priv_verifyRA_supportFunc(linsys,params,options,spec);

    % stl
    case 'stl:kochdumper'
        varargout = cell(1,2);
        [res,varargout{:}] = priv_verifySTL_kochdumper(linsys,params,options,spec);
    case 'stl:seidl'
        varargout = cell(1,0);
        res = verifySTL_seidl(linsys,params,options,spec);

    otherwise
        throw(CORAerror('CORA:noops',linsys,params,options,spec));
end

end


% Auxiliary functions -----------------------------------------------------

function [verifyAlg,options] = aux_selectVerifyAlgorithm(linsys,params,options,spec)

% check if algorithm is given
if isfield(options,'verifyAlg')
    verifyAlg = options.verifyAlg;
    options = rmfield(options,'verifyAlg');
elseif isa(spec,'stl')
    verifyAlg = 'stl:seidl';
elseif isa(spec,'specification')
    % the supportFunc algorithm is definitely faster for specifications
    % given as halfspaces. For others, one would need check which algorithm
    % is faster. Default is 'reachavoid:supportFunc' for now ...
    verifyAlg = 'reachavoid:supportFunc';
else
    throw(CORAerror('CORA:wrongValue','options.verifyAlg',...
        'Verification algorithm not specified and unable to determine one automatically.'))
end


if CHECKS_ENABLED
    % validate chosen verification algorithm
    possibleValues = {
        'reachavoid:zonotope','reachavoid:supportFunc', ...
        'stl:kochdumper','stl:seidl'
    };
    if ~ismember(verifyAlg,possibleValues)
        validrange = ['''', strjoin(possibleValues,"', '"), ''''];
        throw(CORAerror("CORA:wrongValue",'options.verifyAlg',validrange))
    end

    if isa(spec,'stl') && ~startsWith(verifyAlg,'stl:')
        throw(CORAerror("CORA:wrongValue",'options.verifyAlg',...
            'Given stl specifications need an algorithm designed for stl.'))
    end

    % check specification
    if isa(spec,'specification')
        for i=1:length(spec)
            % check each specification
            if strcmp(spec(i).type,'logic')
                if ~startsWith(verifyAlg,'stl:')
                    throw(CORAerror("CORA:wrongValue",'options.verifyAlg', ...
                                        'Given logic specifications need an algorithm designed for stl.'))
                end
                
            else
                % only reachavoid is implemented
                if ~startsWith(verifyAlg,'reachavoid:')
                    throw(CORAerror("CORA:wrongValue",'options.verifyAlg', ...
                        'Given reach-avoid specifications need an algorithm designed for reachavoid.'))
                end
            end
        end
    end
end


end

% ------------------------------ END OF CODE ------------------------------
