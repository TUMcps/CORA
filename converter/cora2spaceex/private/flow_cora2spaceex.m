function flow_cora2spaceex(Obj, location, docNode)
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
%
% Author:        Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---
%
%------------- BEGIN CODE --------------

x = sym('x',[Obj.dim,1]);
u = sym('u',[Obj.nrOfInputs,1]);

% ---------------------------------------------------------------------
%   Linear differential equation
%   eq: x' = A x + B u + c

if isa(Obj,'linearSys')
    
    A   = Obj.A;
    B   = Obj.B;
    
    if ~isempty(Obj.c)
       eq = A*x + B*u + Obj.c; 
    else
       eq = A*x + B*u; 
    end
end

% ---------------------------------------------------------------------
%   Nonlinear differential equation

if isa(Obj,'nonlinearSys')
    
    eq = Obj.mFile(x,u);
    
end

% ---------------------------------------------------------------------
% convert the equation to a character sequence
eqs ='';
for idx = 1:Obj.dim
    x = sprintf('x%d''', idx); % first derivative of x_n
    eq_c = char(eq(idx));
    eq_c = [x, ' == ', eq_c];
    if idx ~= Obj.dim;  eq_c = [eq_c,newline,' & ']; end
    eqs = [eqs,eq_c];
end

%Add the element node (flow), for the parent element (location) and
%set the equation attribute.
flow = docNode.createElement('flow');
flow.appendChild(docNode.createTextNode(eqs));
location.appendChild(flow);

end

%------------- END OF CODE --------------