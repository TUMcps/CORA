function res = test_fourier()
% test_fourier - tests whether fourier can be executed on the current
%    operating system
%
% Syntax:
%    res = test_fourier()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       03-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    % example taken from fourier.m
    H = [-1 -1 -1 0; 3 -1 -1 1; -1 3 -1 2; -1 -1 3 3]; 
    h = fourier(H,[2 3]);

catch ME
    if strcmp(ME.identifier,'MATLAB:TooManyInputs')
        CORAwarning("CORA:global", ...
            ['CORA does not have a compiled file of the Fourier-Motzkin elimination algorithm for your distribution. ' ...
            'Please let us know which distribution you are using and try to install it via the tbxmanager: tbxmanager install fourier.'])
        rethrow(ME)
    end
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
