function res = testLong_polyZonotope_subs
% testLong_polyZonotope_subs - unit test function for composition
%
% Syntax:
%    res = testLong_polyZonotope_subs
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       23-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nTests = 50;
for i=1:nTests
    n = randi(30);
    pZ = polyZonotope.generateRandom('Dimension',n);
    ne_max = size(pZ.E,1);
    % produce random ids for pZ
    id = randperm(ne_max,ne_max)'-ceil(ne_max/2);
    id_subs = id(unique(randi([1,ne_max],ne_max,1),'stable'));
    ns = length(id_subs);
    % generate random pZ to insert
    pZ_in = polyZonotope.generateRandom('Dimension',ns);
    nin = size(pZ_in.E,1);
    id_in = randperm(4*nin,nin)'-2*nin;
    % remove rest matrices
    pZ = polyZonotope(pZ.c,pZ.G,zeros(n,0),pZ.E,id);
    pZ_in = polyZonotope(pZ_in.c,pZ_in.G,zeros(ns,0),pZ_in.E,id_in);
    
    % compute result
    id_rem = setdiff(pZ.id,id_subs,'stable');
    pZ_res = subs(pZ,pZ_in,id_subs);
    % get id indices
    [~,~,ii_rem] = intersect(id_rem,pZ_res.id,'stable');
    [~,~,ii_in] = intersect(id_in,pZ_res.id,'stable'); 
    assert(isempty(setdiff([ii_rem;ii_in],(1:length(pZ_res.id)))),...
            'Error: Remaining and subs id should be all');
    
    % compare with ...
    
    f_pz = fhandle(pZ,{id_rem,id_subs});
    f_in = fhandle(pZ_in,{id_in});
    f_res = @(xr,xin)f_pz(xr,f_in(xin));
    
    % ... by testing random points
    N = 50;
    % could be any real vector really
    B = 2*rand(size(pZ_res.E,1),N)-1;
    for k=1:N
        val_res = resolve(pZ_res,B(:,k));
        fval_res = f_res(B(ii_rem,k),B(ii_in,k));
        a_max = max(abs(val_res),abs(fval_res));
        ind_0 = a_max==0;
        a_max(ind_0) = 1;
        if any(abs(val_res-fval_res)./a_max>1e-6)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
