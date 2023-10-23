function flow_cora2spaceex(obj, location, docNode)
% flow_cora2spaceex - 
%
% Syntax:
%   flow_cora2spaceex(Obj, location, docNode)
%
% Inputs:
%    Obj -
%    location -
%    docNode - 
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

x = sym('x',[obj.dim,1]);
u = sym('u',[obj.nrOfInputs,1]);

% ---------------------------------------------------------------------
%   Linear differential equation
%   eq: x' = A x + B u + c

if isa(obj,'linearSys')
    
    A   = obj.A;
    B   = obj.B;
    
    if ~isempty(obj.c)
       eq = A*x + B*u + obj.c; 
    else
       eq = A*x + B*u; 
    end
end

% ---------------------------------------------------------------------
%   Nonlinear differential equation

if isa(obj,'nonlinearSys')
    
    eq = obj.mFile(x,u);
    
end

% ---------------------------------------------------------------------
% convert the equation to a character sequence
eqs ='';
for idx = 1:obj.dim
    x = sprintf('x%d''', idx); % first derivative of x_n
    eq_c = char(eq(idx));
    eq_c = [x, ' == ', eq_c];
    if idx ~= obj.dim;  eq_c = [eq_c,newline,' & ']; end
    eqs = [eqs,eq_c];
end

%Add the element node (flow), for the parent element (location) and
%set the equation attribute.
flow = docNode.createElement('flow');
flow.appendChild(docNode.createTextNode(eqs));
location.appendChild(flow);

end

% ------------------------------ END OF CODE ------------------------------
