function [res_t1, res_t2] = rescale_dim( t1, t2 )
% rescale_dim creates a new dimention space containing all unique variables
%
% Syntax:  
%    [res_t1, res_t2] = rescale_dim( t1, t2 )
%
% Inputs:
%    t1, t2 - taylm
%
% Outputs:
%    res_t1, res_t2 - rescaled taylm
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      20-April-2016
%               21-July-2016 (DG) the polynomial part is changed to syms
%               18-July-2017 (DG) Multivariable polynomial pack is added
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


x1 = t1.monomials(:, 2:end);
m1 = length( t1.names_of_var );

x2 = t2.monomials(:, 2:end);
m2 = length( t2.names_of_var );

% find unique variables in t1 and t2. They construct a new number of
% dimentions

    
names1 = t1.names_of_var;
names2 = t2.names_of_var;

unique_dim = length(names1);
united_names = {};
counter = 1;

for j = 1:length(names2)
    idx = 0;
    for i = 1:length(names1)
        idx = strcmp(char(names1{i}), char(names2{j})  );
        char(names1{i});
        char(names2{j});
        if idx
            k{j} = i;
            break
        end
    end
    if ~idx
        
        unique_dim = unique_dim + 1;
        united_names{counter} = names2{j};
        counter = counter + 1;
        k{j} = unique_dim;
    end
    
end

% new names of variables
united_names = [names1, flipud(united_names)];

% crate a transformation matrix for t1
T1 = [eye(m1); zeros(unique_dim - m1, m1)];

% crate a transformation matrix for t2
T2 = zeros(unique_dim, m2);

for i = 1:m2
    T2(k{i}, i) = 1;
end

% pack vectors of monomials to ranks
x1 = T1 * x1';
x2 = T2 * x2';

t1.monomials = hashFunction( x1' );
t1.names_of_var = united_names;

t2.monomials = hashFunction( x2' );
t2.names_of_var = united_names;

res_t1 = t1;
res_t2 = t2;

%------------- END OF CODE --------------