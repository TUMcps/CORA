function res = isCanonicalForm(linsys)
% isCanonicalForm - checks if a linear system is formulated in canonical
%    form: x' = Ax + u, y = Cx + v
%
% Syntax:
%    res = isCanonicalForm(linsys)
%
% Inputs:
%    linsys - linearSys object
%
% Outputs:
%    res - true/false whether system is in canonical form
%
% Example:
%    linsys = linearSys(eye(2));
%    res = isCanonicalForm(linsys);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% for readability reasons, perform the checks in batches

% check dimensions
if linsys.nrOfInputs ~= linsys.nrOfDims || linsys.nrOfNoises ~= linsys.nrOfOutputs
    res = false;
    return
end

% ensure that offset vectors, throughput/disturbance matrix are all-zero
if any(linsys.c) || any(linsys.k) || any(linsys.E,'all') || any(linsys.D,'all')
    res = false;
    return
end

% ensure that input/noise matrix are identity
if ~all(linsys.B == eye(linsys.nrOfDims),'all') || ...
        ~all(linsys.C == eye(linsys.nrOfOutputs),'all')
    res = false;
    return
end

% all checks ok -> system is in canonical form
res = true;

% ------------------------------ END OF CODE ------------------------------
