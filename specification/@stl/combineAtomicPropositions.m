function [phi,aps] = combineAtomicPropositions(phi)
% combineAtomicPropositions - combine all atomic propositions of the
%                             formula to geometric set representations
% 
% Syntax:
%    phi = combineAtomicPropositions(phi)
%
% Inputs:
%    phi - STL formula
%
% Outputs:
%    phi - the formula with atomic propositions combined
%    aps - a map of all atomic propositions to their geometric sets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Benedikt Seidl
% Written:       11-May-2023
% Last update:   08-February-2024 (FL, use interval property of stl instead of from and to)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

aps = containers.Map;
phi = aux_inner(phi);


% Auxiliary functions -----------------------------------------------------

function phi = aux_inner(phi)

    if strcmp(phi.type, 'variable') && isa(phi.lhs, 'atomicProposition') % ---
        aps(formattedDisplayText(phi)) = phi.lhs;

    elseif ismember(phi.type,{'<=','<','>','>='}) % ---
        % combine complete subformula to geometric set
        ap = atomicProposition(convert2set(phi));
        aps(formattedDisplayText(phi)) = ap;
        phi = stl(char(formattedDisplayText(phi)), ap);

    elseif strcmp(phi.type, 'until') % ---
        lhs = aux_inner(phi.lhs);
        rhs = aux_inner(phi.rhs);
        phi = until(lhs, rhs, phi.interval);

    elseif strcmp(phi.type, 'release') % ---
        lhs = aux_inner(phi.lhs);
        rhs = aux_inner(phi.rhs);
        phi = release(lhs, rhs, phi.interval);

    elseif strcmp(phi.type, 'finally') % ---
        lhs = aux_inner(phi.lhs);
        phi = finally(lhs, phi.interval);

    elseif strcmp(phi.type, 'globally') % ---
        lhs = aux_inner(phi.lhs);
        phi = globally(lhs, phi.interval);

    elseif strcmp(phi.type, 'next') % ---
        lhs = aux_inner(phi.lhs);
        phi = next(lhs, phi.from);
        
    elseif strcmp(phi.type, '&') % ---
        lhs = aux_inner(phi.lhs);
        rhs = aux_inner(phi.rhs);
        phi = lhs & rhs;

    elseif strcmp(phi.type, '|') % ---
        lhs = aux_inner(phi.lhs);
        rhs = aux_inner(phi.rhs);
        phi = lhs | rhs;

    elseif strcmp(phi.type, '~') % ---
        phi = ~ aux_inner(phi.lhs);
    end

end

end

% ------------------------------ END OF CODE ------------------------------
