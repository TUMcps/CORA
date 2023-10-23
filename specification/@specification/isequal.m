function res = isequal(spec1,spec2,varargin)
% isequal - checks if two specification objects are equal
%
% Syntax:
%    res = isequal(spec1,spec2)
%    res = isequal(spec1,spec2,tol)
%
% Inputs:
%    spec1 - specification object
%    spec2 - specification object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    set = halfspace([1 0],2);
%    spec1 = specification(set,'unsafeSet');
%    spec2 = specification(set,'safeSet');
%    res = isequal(spec1,spec2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       29-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default value
tol = setDefaultValues({1e-12},varargin);

% check input arguments
inputArgsCheck({{spec1,'att','specification'},...
                {spec2,'att','specification'},...
                {tol,'att','numeric',{'nonnegative,','scalar'}}});

% specification objects have to have same length
nrSpecs = length(spec1);
if nrSpecs ~= length(spec2)
    res = false; return
end

% loop over individual specifications trying to find matches
spec2Matched = false(nrSpecs,1);
for i=1:nrSpecs
    for j=1:nrSpecs
        if ~spec2Matched(j) && aux_sameSpecs(spec1(i),spec2(j),tol)
            % j-th specification of spec2 = i-th specification of spec1
            spec2Matched(j) = true; break
        end
    end
    if nnz(spec2Matched) < i
        % could not find a match for i-th specification in spec1
        res = false; return
    end
end

% all checks ok
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_sameSpecs(spec1,spec2,tol)

% assume true
res = true;

% type has to be the same
if ~strcmp(spec1.type,spec2.type)
    res = false; return
end

% time interval where specifications are active have to be the same
if ~isequal(spec1.time,spec2.time,tol)
    res = false; return
end

% either both specs have a defined location or none
if xor(isempty(spec1.location),isempty(spec2.location))
    res = false; return
end

if ~isempty(spec1.location)
    % locations have to have equal size
    if ~all(size(spec1.location) == size(spec2.location))
        res = false; return
    end

    if ~iscell(spec1.location)
        % location for hybrid automata
        if ~all(spec1.location == spec2.location)
            res = false; return
        end
    else
        % location for parallel hybrid automata (loop over cells)
        for i=1:length(spec1.location)
            if ~all(spec1.location{i},spec2.location{i})
                res = false; return
            end
        end
    end
end

% safe sets/unsafe sets/invariants: compare sets (same type ensured above)
if strcmp(spec1.type,'safeSet') || strcmp(spec1.type,'unsafeSet') ...
        || strcmp(spec1.type,'invariant')
    if ~isequal(spec1.set,spec2.set,tol)
        res = false; return
    end
end

% --- currently cannot compare stl and costum

% compare stl formula
if strcmp(spec1.type,'logic') %&& ~isequal(spec1.set,spec2.set,tol)
    throw(CORAerror('CORA:notSupported','isequal not supported for stl formulae.'));
end

% compare costum (function handle)
if strcmp(spec1.type,'costum')
    throw(CORAerror('CORA:notSupported','isequal not supported for custom functions.'));
end

end

% ------------------------------ END OF CODE ------------------------------
