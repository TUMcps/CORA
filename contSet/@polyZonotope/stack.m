function pZ = stack(varargin)
% stack - stacks polyZonotopes to retain dependencies
%
% Syntax:  
%    PZ = stack(pZ1,pZ2,...)
%    PZ = subs(pZ,[],...)
%
% Inputs:
%    pZi - polyZonotope object or []
%
% Outputs:
%    PZ - stacked polyZonotope object
%
% Example: 
%    pZ1 = polyZonotope(0,[1,1],[],[1,0;2,1],[1;2]);
%    pZ2 = polyZonotope([0;1],eye(2),zeros(2,0),[2,0;0,1],[1;3]);
%    res = stack(pZ1,pZ2);
%    print(res,{[1;2],3},{'b','r'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cartProd

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    pZ = [];
    i = 1;
    while i<=nargin
        inp1 = varargin{i};
        inp2 = [];
        if i+1<=nargin
            inp2 = varargin{i+1};
        end
        if (~isa(inp1,'polyZonotope') && ~isempty(inp1)) || ...
           (~isa(inp2,'polyZonotope') && ~isempty(inp2))
            error('Wrong input arguments');
        end
        if ~isempty(inp1) || ~isempty(inp2)
            pZ = stack2(pZ,stack2(inp1,inp2));
        end
        i = i+2;
    end
end
%%------%%
function pZ = stack2(pZ1,pZ2)
    if isempty(pZ1)
        pZ = pZ2;
        return;
    elseif isempty(pZ2)
        pZ = pZ1;
        return;
    end
    n1 = size(pZ1.c,1);
    n2 = size(pZ2.c,1);
    pZ1t = polyZonotope([pZ1.c;zeros(n2,1)],...
                        [pZ1.G;zeros(n2,size(pZ1.G,2))],...
                        [pZ1.Grest;zeros(n2,size(pZ1.Grest,2))],...
                        pZ1.expMat,pZ1.id);
    pZ2t = polyZonotope([zeros(n1,1);pZ2.c],...
                        [zeros(n1,size(pZ2.G,2));pZ2.G],...
                        [zeros(n1,size(pZ2.Grest,2));pZ2.Grest],...
                        pZ2.expMat,pZ2.id);
    % check if one of them is only a vertex
    if isempty(pZ1t.G)
        pZ = polyZonotope(pZ2t.c+pZ1t.c,pZ2t.G,pZ2t.Grest,pZ2t.expMat,pZ2t.id);
    elseif isempty(pZ2t.G)
        pZ = polyZonotope(pZ1t.c+pZ2t.c,pZ1t.G,pZ1t.Grest,pZ1t.expMat,pZ1t.id);
    else
        pZ = exactPlus(pZ1t,pZ2t);
    end
end
%------------- END OF CODE --------------