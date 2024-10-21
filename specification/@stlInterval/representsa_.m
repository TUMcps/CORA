function [res,S] = representsa_(I,type,tol,varargin)
% representsa_ - checks if an STL interval can also be represented
%    by a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(I,type,tol)
%    [res,S] = representsa_(I,type,tol)
%
% Inputs:
%    I - stlInterval object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Florian Lercher
% Written:       09-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(I,type);
else
    [empty,res,S] = representsa_emptyObject(I,type);
end
if empty; return; end

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'stlInterval'
        % obviously true
        res = true;
        if nargout == 2
            S = I;
        end
    otherwise
        if I.leftClosed && I.rightClosed
            closedInterval = interval(I);
            if nargout == 1
                res = representsa(closedInterval,type,tol,varargin{:});
            else
                [res,S] = representsa(closedInterval,type,tol,varargin{:});
            end
        else
            % all other set representations are closed sets,
            % so they cannot be represented by a (half-)open interval
            res = false;
        end
end

% ------------------------------ END OF CODE ------------------------------
