function E = ellipsoid(I, varargin)
% ellipsoid - Converts an interval object into an ellipsoid object
%
% Syntax:
%    E = ellipsoid(I)
%    E = ellipsoid(I,mode)
%
% Inputs:
%    I - interval object
%    mode - str, 'outer' or 'inner'
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    I = interval([1;-1], [2; 1]);
%    E = ellipsoid(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Victor Gassmann, Tobias Ladner
% Written:       15-October-2019 
% Last update:   18-April-2024 (TL, bug fix, almost degenerate dimensions)
%                18-April-2024 (TL, add 'inner' case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. preprocess
[I,mode] = aux_preprocess(I,varargin);

% 2. compute ellipsoid
if startsWith(mode, 'outer')
    % outer approximation

    % avoid numerical issues, enlarge very small non-degenerate dimensions
    % TL, seems to occur below radius=5e-3 (squared in zonotope/ellipsoid)

    radMin = 5e-3;
    radCur = rad(I);

    factor = 1.*(radCur>=radMin) + (radMin./radCur).*(radCur<radMin);
    factor(isinf(factor)) = 1;
    I = enlarge(I,factor);

    % convert interval to zonotope (exact computation very efficient because
    % generator matrix is square)
    E = ellipsoid(zonotope(I),'outer:exact');

elseif startsWith(mode, 'inner')
    % inner approximation

    % convert interval to zonotope (exact computation very efficient because
    % generator matrix is square)
    E = ellipsoid(zonotope(I),'inner:exact');
end

end


% Auxiliary functions -----------------------------------------------------

function [I,mode] = aux_preprocess(I,givenvalues)
    
    % set default values
    [mode] = setDefaultValues({'outer'},givenvalues);

    % check input args
    inputArgsCheck({ ...
       {I,'att','interval'}; ...
       {mode,'str',{'outer','inner'}}; ...
    })

end

% ------------------------------ END OF CODE ------------------------------
