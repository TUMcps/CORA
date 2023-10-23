function pZ = stack(varargin)
% stack - stacks polynomial zonotopes to retain dependencies (in some sense
%           a cartesian product which retains dependencies)
%
% Syntax:
%    pZ = stack(pZ1,pZ2,...)
%
% Inputs:
%    pZi1,pZ2,... - polyZonotope object or []
%
% Outputs:
%    pZ - stacked polyZonotope object
%
% Example: 
%    pZ1 = polyZonotope(0,[1,1],[],[1,0;2,1],[1;2]);
%    pZ2 = polyZonotope([0;1],eye(2),zeros(2,0),[2,0;0,1],[1;3]);
%    pZ = stack(pZ1,pZ2);
%    print(pZ,'Ids',{[1;2],3},'Vars',{'b','r'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cartProd

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   04-July-2022 (VG, input checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pZ = [];
% add [] at the end if not even
if mod(nargin,2)>0
    varargin{end+1} = [];
end
for i=1:2:nargin
    % extract and check inputs (if not empty)
    inp1 = varargin{i};
    if ~isempty(inp1)
        inputArgsCheck({{inp1,'att','polyZonotope'}});
    end
    inp2 = varargin{i+1};
    if ~isempty(inp2)
        inputArgsCheck({{inp2,'att','polyZonotope'}});
    end

    if ~isempty(inp1) || ~isempty(inp2)
        % recursively build up pZ
        pZ = aux_stack2(pZ,aux_stack2(inp1,inp2));
    end
end
end


% Auxiliary functions -----------------------------------------------------

function pZ = aux_stack2(pZ1,pZ2)
    % if either is empty, simply return the other
    if representsa_(pZ1,'emptySet',eps)
        pZ = pZ2;
        return;
    elseif representsa_(pZ2,'emptySet',eps)
        pZ = pZ1;
        return;
    end
    
    % extract dimensions
    n1 = dim(pZ1);
    n2 = dim(pZ2);
    % extend both to final dimensions (n1+n2)
    pZ1t = polyZonotope([pZ1.c;zeros(n2,1)],...
                        [pZ1.G;zeros(n2,size(pZ1.G,2))],...
                        [pZ1.GI;zeros(n2,size(pZ1.GI,2))],...
                        pZ1.E,pZ1.id);
    pZ2t = polyZonotope([zeros(n1,1);pZ2.c],...
                        [zeros(n1,size(pZ2.G,2));pZ2.G],...
                        [zeros(n1,size(pZ2.GI,2));pZ2.GI],...
                        pZ2.E,pZ2.id);

    % check if one of them is only a vertex
    if isempty(pZ1t.G)
        pZ = polyZonotope(pZ2t.c+pZ1t.c,pZ2t.G,pZ2t.GI,pZ2t.E,pZ2t.id);
    elseif isempty(pZ2t.G)
        pZ = polyZonotope(pZ1t.c+pZ2t.c,pZ1t.G,pZ1t.GI,pZ1t.E,pZ1t.id);
    else
        pZ = exactPlus(pZ1t,pZ2t);
    end

end

% ------------------------------ END OF CODE ------------------------------
