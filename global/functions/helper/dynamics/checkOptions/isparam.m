function res = isparam(name)
% isparam - determines whether a given name is part of the model parameters
%    or the algorithm parameters
%
% Syntax:
%    res = isparam(name)
%
% Inputs:
%    name - field name of params struct or options struct
%
% Outputs:
%    res - true if param, false if option
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: splitIntoParamsOptions

% Authors:       Mark Wetzlinger
% Written:       21-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

switch name
    case {'tStart','tFinal','R0','U','u','W','V','startLoc','finalLoc',...
            'inputCompMap','x0','paramInt','y0guess','tu','safeSet','unsafeSet'}
        % note: move safeSet and unsafeSet to spec!
        res = true;
    otherwise
        % lazy approach: all others are options
        res = false;
%     otherwise
%         throw(CORAerror('CORA:specialError','Unknown field name.'));
end

% ------------------------------ END OF CODE ------------------------------
