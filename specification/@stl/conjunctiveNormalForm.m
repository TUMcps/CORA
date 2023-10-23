function res = conjunctiveNormalForm(obj,varargin)
% conjunctiveNormalForm - convert STL formula to conjunctive normal form
%
% Syntax:
%    res = conjunctiveNormalForm(obj)
%    res = conjunctiveNormalForm(obj,nnf)
%
% Inputs:
%    obj - logic formula (class stl)
%    neg - flag specifying weather the formula is in negation normal form 
%
% Outputs:
%    res - resulting stl formula in conjunctive normal form (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = ~(x(1) < 5 | x(2) < 3) | x(2) > 5;
%    eq_ = conjunctiveNormalForm(eq)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    nnf = false;

    if nargin > 1 
        nnf = varargin{1};
    end

    % convert to negation normal form
    if ~nnf
        obj = negationNormalForm(obj);
    end

    % convert to conjunctive normal form
    while true
        [obj,cnt] = aux_rewrite(obj);
        if cnt == 0
            break;
        end
    end

    res = obj;
end


% Auxiliary functions -----------------------------------------------------

function [res,cnt] = aux_rewrite(obj)
% recursive function to convert into conjunctive normal form by
% distributing & over |

    if strcmp(obj.type,'|')
        
        [lhs,cnt1] = aux_rewrite(obj.lhs);
        [rhs,cnt2] = aux_rewrite(obj.rhs);

        cnt = cnt1 + cnt2;

        if strcmp(lhs.type,'&') && strcmp(rhs.type,'&')
            res = (lhs.lhs | rhs.lhs) & (lhs.lhs | rhs.rhs) & ...
                            (lhs.rhs | rhs.lhs) & (lhs.rhs | rhs.rhs);
            cnt = cnt + 1;
        elseif strcmp(lhs.type,'&')
            res = (rhs | lhs.lhs) & (rhs | lhs.rhs);
            cnt = cnt + 1;
        elseif strcmp(obj.rhs.type,'&')
            res = (lhs | rhs.lhs) & (lhs | rhs.rhs);
            cnt = cnt + 1;
        else
            res = lhs | rhs;
        end
    elseif strcmp(obj.type,'&')

        [lhs,cnt1] = aux_rewrite(obj.lhs);
        [rhs,cnt2] = aux_rewrite(obj.rhs);

        cnt = cnt1 + cnt2;

        res = lhs & rhs;
    else
        res = obj;
        cnt = 0;
    end
end

% ------------------------------ END OF CODE ------------------------------
