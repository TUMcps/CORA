function assignment_cora2spaceex(tran, docNode, reset)
% assignment_cora2spaceex - 
%
% Syntax:
%    assignment_cora2spaceex(tran, docNode, reset)
%
% Inputs:
%    tran -
%    docNode - 
%    reset - 
%
% Outputs:
%    -
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty reset function (0x0 transition)
if isempty(fields(reset))
    eqs = '';
    
elseif isfield(reset,'A')
    % linear reset

    A = reset.A;
    b = reset.c;
    
    dim = size(A,1);
    
    x = sym('x',[dim,1]);
    
    eq = A*x + b;
    
    % convert the equation to a character sequence
    eqs ='';
    for idx = 1:dim
        x = sprintf('x%d''', idx); % first derivative of x_n
        eq_c = char(eq(idx));
        eq_c = [x, ' := ', eq_c];
        if idx > 1;  eq_c = [newline,' & ', eq_c]; end
        eqs = [eqs,eq_c];
    end
elseif isfield(reset,'f')
    throw(CORAerror('CORA:notSupported',...
        'cora2spaceex does not support nonlinear reset functions.'));
end

%Add the element node (assignment), for the parent element (transition) and
%set the equation attribute.
assignment = docNode.createElement('assignment');
assignment.appendChild(docNode.createTextNode(eqs));
tran.appendChild(assignment);

end

% ------------------------------ END OF CODE ------------------------------
