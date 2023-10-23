function res = test_polyZonotope_cubMap
% test_polyZonotope_cubMap - unit test function for the cubic 
%    multiplication of polynomial zonotopes with a tensor
%
% Syntax:
%    res = test_polyZonotope_cubMap
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

% Authors:       Niklas Kochdumper
% Written:       17-August-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TEST cubic multiplication

% define polynomial zonotope
pZ = polyZonotope([0;1],[1 -1;2 0],[],eye(2));

% define third-order tensor
temp = [1 -1; 0 2];
T{1,1} = temp;
T{1,2} = temp;
T{2,1} = temp;
T{2,2} = temp;

% compute cubic map
pZres = cubMap(pZ,T);

% define ground truth
temp = [2 13 -1 28 -4 21 -7 3 -1];
Z_ = [temp;temp];
c_ = Z_(:,1);
G_ = Z_(:,2:end);
E_ = [1 0 2 1 3 2 1 0;...
           0 1 0 1 0 1 2 3];
       
% check for correctness
if ~all(withinTol(pZres.c,c_)) || any(size(pZres.G)-size(G_)) ...
        || any(size(pZres.E)-size(E_))
    throw(CORAerror('CORA:testFailed'));
else
    for i = 1:size(E_,2)    

        ind = ismember(pZres.E',E_(:,i)','rows');  
        ind_ = find(ind > 0);

        if isempty(ind_)
            throw(CORAerror('CORA:testFailed'));
        elseif ~all(pZres.G(:,ind_(1)) == G_(:,i))
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% second test with extreme point

% pZ
c = [-0.1383; -0.2355];
G = [; ...
    0.4654, -0.1336, 0.2960, 0.3770, 0.2401, -0.1173, -0.1465, 0.0956; ...
    -0.4872, 0.3647, -0.3217, -0.1125, 0.2805, -0.3750, 0.1550, 0.3522];
GI = [-0.3132, -0.1837, -0.1806; ...
    0.4952, -0.3656, -0.4542];
E = [; ...
    1, 0, 0, 2, 2, 8, 6, 3; ...
    0, 1, 0, 7, 1, 3, 4, 2; ...
    0, 0, 1, 3, 2, 1, 4, 9];
id = [1; 2; 3];
pZ = polyZonotope(c, G, GI, E, id);

% extreme point
p0 = [-1.0527; 3.5282];

% third-order tensor
T = {; ...
    [-0.3274, -0.3096; 0.0991, -0.3295], [0.1074, 0.0649; -0.4000, 0.3806]; ...
    [-0.4289, 0.1701; 0.3544, -0.4907], [0.0067, -0.0605; 0.3368, -0.2607]; ...
    };

% compute cubic map
pZ_res = cubMap(pZ, T);
p0_res = cubMapPoint(p0, p0, p0, T);

% figure; 
% % input
% subplot(1, 2, 1); hold on;
% plot(pZ,[1 2],'Splits',14);
% scatter(p0(1,:),p0(2,:),'.k')
% % output
% subplot(1, 2, 2); hold on;
% plot(pZ_res,[1 2],'Splits',14);
% scatter(p0_res(1,:),p0_res(2,:),'.k')

if ~containsPointSet(pZ_res, p0_res, [], 30)
    throw(CORAerror('CORA:testFailed'));
end


res = true;

% ------------------------------ END OF CODE ------------------------------
