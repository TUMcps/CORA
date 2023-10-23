function phi = desugar(phi)
% desugar - remove syntactic sugar from the formula
%
% The only operators appearing in the resulting formula are until,
% conjunction, disjunction and negation.
%
% Syntax:
%    phi = desugar(phi);
%
% Inputs:
%    phi - STL formula
%
% Outputs:
%    phi - the desugared formula
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Benedikt Seidl
% Written:       11-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

arguments
    phi stl
end

if strcmp(phi.type, 'until')
    phi = until(desugar(phi.lhs), desugar(phi.rhs), interval(phi.from, phi.to));
elseif strcmp(phi.type, 'release')
    phi = desugar(~until(~phi.lhs, ~phi.rhs, interval(phi.from, phi.to)));
elseif strcmp(phi.type, 'finally')
    phi = desugar(until(stl(true), phi.lhs, interval(phi.from, phi.to)));
elseif strcmp(phi.type, 'globally')
    phi = desugar(~finally(~phi.lhs, interval(phi.from, phi.to)));
elseif strcmp(phi.type, 'next')
    phi = desugar(finally(phi.lhs, interval(phi.from, phi.from)));
elseif strcmp(phi.type, '&')
    phi = desugar(phi.lhs) & desugar(phi.rhs);
elseif strcmp(phi.type, '|')
    phi = desugar(phi.lhs) | desugar(phi.rhs);
elseif strcmp(phi.type, '~')
    phi = ~ desugar(phi.lhs);
end

end

% ------------------------------ END OF CODE ------------------------------
