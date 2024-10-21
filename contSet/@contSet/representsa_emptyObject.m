function [empty,res,S_conv] = representsa_emptyObject(S,type)
% representsa_emptyObject - checks if a contSet class object is fully empty
%    and converts the set to an instance of class 'type'
%
% Syntax:
%    [empty,res] = representsa_emptyObject(S,type)
%    [empty,res,S_conv] = representsa_emptyObject(S,type)
%
% Inputs:
%    S - contSet object
%    type - contSet class
%
% Outputs:
%    empty - true/false whether S is fully empty
%    res - true/false whether S can be represented by 'type'
%    S_conv - object of class 'type', converted from S
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   01-January-2024 (MW, update fully empty polytopes)
%                14-July-2024 (MW, support polytopes in V representation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if set is a fully empty object
empty = isemptyobject(S);
% default: no converted set
res = [];
S_conv = [];

% dimension
n = dim(S);

% only further information is S is fully empty
if empty

    % self-checking
    if isa(S,type)
        res = true;
        if nargout == 3 
            S_conv = eval([type, '.empty(', num2str(n), ')']);
        end
        return
    end

    switch type
        case {'origin','point','parallelotope'}
            % fully empty objects cannot be any single point or a parallelotope
            res = false;
            
        case 'fullspace'
            % fully empty polytopes in halfspace representation represent
            % R^n (=fullspace); ensure that V representation is not given
            res = isa(S,'polytope') && ...
                ( (S.isHRep.val && isempty(S.b) && isempty(S.be)) ...
                || (S.isVRep.val && ~isempty(S.V) && all(all(isinf(S.V)))) );
            if nargout == 3 && res
                S_conv = fullspace(dim(S));
            end
    
        case 'hyperplane'
            res = true;

        case 'interval'
            res = true;
            if nargout == 3
                if isa(S,'polytope') && ...
                ( (S.isHRep.val && isempty(S.b) && isempty(S.be)) ...
                || (S.isVRep.val && ~isempty(S.V) && all(all(isinf(S.V)))) )
                    % fullspace
                    S_conv = interval(-Inf(dim(S),1),Inf(dim(S),1));
                else
                    % empty set
                    S_conv = interval(zeros(dim(S),0));
                end
            end

        case 'polytope'
            res = true;
            if nargout == 3
                % all other fully empty set reps represent the empty set
                S_conv = polytope.empty(dim(S));
            end
    
        otherwise
            % all fully empty objects represent the empty set (except for

            % polytopes and spectrahedral shadow); all sets can represent
            % the empty set
            res = dim(S) == 0 || (~isa(S,'polytope') && ~isa(S,'spectraShadow')) || ...
                (isa(S,'polytope') && S.isVRep.val && isempty(S.V));

            if nargout == 3 && res
                S_conv = eval([type, '.empty(', num2str(n), ')']);
            end
    end
end

% ------------------------------ END OF CODE ------------------------------
