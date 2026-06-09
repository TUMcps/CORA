function [specSet,phi,u] = preprocessFalsify(params,spec)
% preprocessFalsify - preprocess the paramters and specifications for a 
%    falsificatoin problem 
%
% Syntax:
%    [specSet,phi,u,N] = preprocessFalsify(params,spec)
%
% Inputs:
%    params - parameter defining the reachability problem
%    spec - object of class specification (reach-avoid) or stl
%
% Outputs:
%    specSet - reach-avoid specification (class specification)
%    phi - temporal logic formula (class stl)
%    u - time-varying input signal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       23-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % split into temporal logic and other specifications
    [spec,specLogic] = splitLogic(spec);

    % extract a list of unsafe and safe sets from the specification object
    specSet = [];

    for i = 1:length(spec)
        if ismember(spec(i).type,{'unsafeSet','safeSet'})

            % convert set to a polytope with unit halfspace normal vectors
            S = polytope(spec(i).set);

            if ~S.isHRep.val
                constraints(S);
            end

            S = normalizeConstraints(S,'A');

            % add the time interval where the specification is active
            if ~representsa(spec(i).time,'emptySet')
                t = spec(i).time;
            else
                t = interval(params.tStart,params.tFinal);
            end

            specSet = [specSet;specification(S,spec(i).type,t)];
        else
            throwAsCaller(CORAerror('CORA:specialError',...
           ['Specifications of type ',spec(i).type,' are not supported']));
        end
    end

    % extract temporal logic formula
    phi = [];

    for i = 1:length(specLogic)
        phi = phi & specLogic(i).set;
    end

    if ~isempty(phi) && params.tStart ~= 0
        phi = combineNext(next(phi,-params.tStart));
    end

    if ~isempty(phi) && ~isfield(params,'u')
        params.tFinal = max(params.tFinal-params.tStart,maximumTime(phi));
    end

    if ~isempty(phi) && maximumTime(phi) > params.tFinal + eps
        throwAsCaller(CORAerror('CORA:specialError',...
                'Maximum time for STL formula exceeds params.tFinal'));
    end

    % determine initial number of time steps
    if ~isfield(params,'u')
        if isfield(params,'uTransVec')
            u = params.uTransVec;
        else
            u = params.uTrans;
        end
    else
        u = params.u;
    end

% ------------------------------ END OF CODE ------------------------------
