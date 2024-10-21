function SpS_out = mtimes(factor1, factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix
%    with a spectrahedral shadow
%
% Syntax:
%    SpS_out = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or scalar
%    factor2 - spectraShadow object 
%
% Outputs:
%    SpS_out - spectrahedral shadow after multiplication with a matrix or
%       scalar
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    M = [2 3;-4 0];
%    MSpS = M * SpS;
%
% Other m-files required: spectraShadow.m
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    %Find a spectraShadow object
    [SpS,M] = findClassArg(factor1,factor2,'spectraShadow');

    if ~isnumeric(M)
        throw(CORAerror('CORA:noops',M,SpS));
    end
    % So M is numeric; check if it's the first or second argument
    if isnumeric(factor1) || (isnumeric(factor2) && isscalar(M))
        % construct spectrahedral shadow
        SpS_out = spectraShadow(SpS.A,M*SpS.c,M*SpS.G);
        
        % Deducing additional properties; computing the rank of M might be
        % expensive, but overall this greatly simplifies the usability of
        % spectrahedral shadows
        SpS_out.emptySet.val = SpS.emptySet.val;
        
        if ~isempty(SpS.center.val)
            % Multiply center by by M
            SpS_out.center.val = M * SpS.center.val;
        end
        
        if isscalar(M) || rank(full(M)) == dim(SpS)
            if ~isempty(SpS.bounded.val) && SpS.bounded.val
                % SpS bounded and M full rank -> SpS_out bounded
                SpS_out.bounded.val = true;
            end
            if ~isempty(SpS.fullDim.val) && SpS.fullDim.val
                % SpS bounded and M full rank -> SpS_out bounded
                SpS_out.fullDim.val = true;
            end
        end
        
    else
        throw(CORAerror('CORA:noops',M,SpS));
    end
    
end

% ------------------------------ END OF CODE ------------------------------
