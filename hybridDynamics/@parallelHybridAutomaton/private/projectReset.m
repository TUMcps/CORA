function res = projectReset(reset,stateBind,dims)
% projectReset - Compute the projection of a reset function into a higher 
%                dimension
%
% Syntax:
%    res = projectReset(reset,stateBind,dims)

% Input:
%    reset - linear map: x -> reset.A*x + reset.b
%    stateBind - binding of the reset's dimensions to overall dimensions of
%                the parallel hybrid automaton
%    dims - overall dimensions of the parallel hybrid automaton
%
% Outputs:
%    res - resulting reset function
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Johann Schöpfer
% Written:      14-June-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % get fields from struct
    A = reset.A;
    b = reset.b;
    
    % default is identity <=> other components remain unaffected
    Aproj = eye(dims,dims);
    bProj = zeros(dims,1);
    
    % reset matrix A
    Aproj(stateBind,stateBind) = A;
    
    % transition vector b
    bProj(stateBind) = b;

    % construct resuling reset object
    res.A = Aproj;
    res.b = bProj;

end

%------------- END OF CODE --------------