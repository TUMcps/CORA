function priv_assignment_cora2spaceex(tran, docNode, reset)
% priv_assignment_cora2spaceex - 
%
% Syntax:
%    priv_assignment_cora2spaceex(tran, docNode, reset)
%
% Inputs:
%    tran - transition object
%    docNode - 
%    reset - linearReset/nonlinearReset object
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

if isa(reset,'linearReset')
    % linear reset

    A = reset.A;
    B = reset.B;
    c = reset.c;
    
    % assumes preStateDim == postStateDim
    n = reset.postStateDim;
    m = reset.inputDim;
    
    x = sym('x',[n,1]);
    u = sym('u',[m,1]);
    
    eq = A*x + c; % + B*u
    
    % convert the equation to a character sequence
    eqs = '';
    for idx = 1:n
        lhs = sprintf('x%d', idx); % first derivative of x_n
        eq_c = char(eq(idx));
        eq_c = [lhs, ' := ', eq_c];
        if idx > 1
            eq_c = [newline,' & ', eq_c];
        end
        eqs = [eqs, eq_c];
    end
elseif isa(reset,'nonlinearReset')
    throw(CORAerror('CORA:notSupported',...
        'cora2spaceex does not support nonlinear reset functions.'));
end

%Add the element node (assignment), for the parent element (transition) and
%set the equation attribute.
assignment = docNode.createElement('assignment');
assignment.appendChild(docNode.createTextNode(eqs));
tran.appendChild(assignment);

% ------------------------------ END OF CODE ------------------------------
