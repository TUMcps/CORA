function res = test_polyZonotope_partZonotope
% test_polyZonotope_partZonotope - unit test function for computing a
%    partial zonotope overapproximation
%
% Syntax:
%    res = test_polyZonotope_partZonotope
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
% Written:       24-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simple test for given polyZonotope
G = [-8,8,-10,3;6,-7,3,-4];
c = [-1;2];
E = [2,1,3,3;0,2,1,2;4,2,2,1];
pZ = polyZonotope(c,G,zeros(2,0),E);
id = pZ.id;
% overapproximate in second id entry
[pZ_res,pZ_gi_c_] = partZonotope(pZ,pZ.id(2));

% remaining id indices
ii = [1,3];
% correct result generators
ii_center = 1;
ii_even = [2,4];
ii_single = 3;
% center first (b_2 = 0) + all even exponents
G1 = [G(:,ii_center),1/2*G(:,ii_even)];
eM1 = E(ii,[ii_center,ii_even]);
pZ_g1 = polyZonotope(c,G1,zeros(2,0),eM1,id(ii));

% all with b_2=1
G2 = G(:,ii_single);
eM2 = E(ii,ii_single);
pZ_g2 = polyZonotope(zeros(2,1),G2,zeros(2,0),eM2,id(ii));

% all with b_2=2
G3 = 1/2*G(:,ii_even);
eM3 = E(ii,ii_even);
pZ_g3 = polyZonotope(zeros(2,1),G3,zeros(2,0),eM3,id(ii));
pZ_gi_c = {pZ_g1,pZ_g2,pZ_g3};
pZ_gi_c_save = pZ_gi_c;

% check if resulting generators polynomial zonotopes are equal
while ~isempty(pZ_gi_c) && ~isempty(pZ_gi_c_)
    
    res_tmp = false;
    for i=1:length(pZ_gi_c_)
        res_tmp = isequal(pZ_gi_c{1},pZ_gi_c_{i});
        if res_tmp
            break;
        end
    end
    res_g = res_tmp;
    pZ_gi_c(1) = [];
    pZ_gi_c_(i) = [];
    if ~res_g
        break;
    end
end
res = res_g && isempty(pZ_gi_c) && isempty(pZ_gi_c_);

% check if pZ_res is correct
id_new = pZ_res.id(~ismember(pZ_res.id,id));
res = res && length([0;id_new])==length(pZ_gi_c_save);
if res
    b = zeros(length(id_new),1);
    pZ_gi_center = resolve(pZ_res,b,id_new);
    if ~isequal(pZ_gi_center,pZ_gi_c_save{1})
        res = false;
    end
    pZ_res = exactPlus(pZ_res,-1*pZ_gi_center);
    pZ_gi_c_save(1) = [];
    if res
        for i=1:length(id_new)
            ind = ismember(id_new,id_new(i));
            b = zeros(length(id_new),1);
            b(ind) = 1;
            pZ_gi_cand = resolve(pZ_res,b,id_new);
            % check if candidate is contained in 
            for j=1:length(pZ_gi_c_save)
                res_tmp = isequal(pZ_gi_cand,pZ_gi_c_save{j});
                if res_tmp
                    break;
                end
            end
            pZ_gi_c_save(j) = [];
            if ~res_tmp
                res = false;
                break;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
