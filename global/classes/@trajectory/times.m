function traj = times(factor1,factor2)
% times - Overloaded '.*' operator for the states of the simulated trajectories
%
% Syntax:
%    traj = times(A,traj)
%
% Inputs:
%    A - numeric vector
%    traj - trajectory object
%
% Outputs:
%    traj - transformed trajectory object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner, Laura Luetzow
% Written:       08-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[traj,A] = findClassArg(factor1, factor2, 'trajectory');

% parse input
if ~isnumeric(A)
    throw(CORAerror('CORA:noops', traj, A))
elseif ~isvector(A) && ~isempty(A)
    throw(CORAerror('CORA:notSupported', 'Multiplied vector has to be a column vector.'))
end

% transform column vector to diagonal matrix and use mtimes
A = diag(A);
traj = A * traj;

end

% ------------------------------ END OF CODE ------------------------------
