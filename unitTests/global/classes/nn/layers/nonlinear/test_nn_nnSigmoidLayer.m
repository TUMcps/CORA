function res = test_nn_nnSigmoidLayer()
% test_nn_nnSigmoidLayer - tests the sigmoid layer
%
% Syntax:
%    res = test_nn_nnSigmoidLayer()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       01-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% init layer
layer = nnSigmoidLayer();

% test values
resvec(end+1) = layer.f(0) == 0.5;
resvec(end+1) = layer.f(inf) == 1;
resvec(end+1) = layer.f(-inf) == 0;

% test name
customName = 'MyLayer';
layer = nnSigmoidLayer(customName);
resvec(end+1) = strcmp(layer.name,customName);

% test region polynomials
reg_polys = layer.reg_polys;
for i=1:length(reg_polys)
    regi = reg_polys(i);

    l = regi.l;
    u = regi.u;

    if isinf(l)
        l = -1000;
    end
    if isinf(u)
        u = 1000;
    end

    % sample points
    xs = linspace(l,u,100);
    ys = layer.f(xs);
    ys_poly = polyval(regi.p,xs);

    % test if error is correct
    resvec(end+1) = all(abs(ys-ys_poly) <= regi.d + eps);


    % test if entire domain is covered
    if i == 1
        % start with -inf
        resvec(end+1) = isequal(regi.l,-inf);
    end
    if i < length(reg_polys)
        % regions connect
        resvec(end+1) = isequal(regi.u,reg_polys(i+1).l);
    end
    if i == length(reg_polys)
        % end with +inf
        resvec(end+1) = isequal(regi.u,+inf);
    end
end

% check evaluate

% check point
x = [1;2;3;4];
y = layer.evaluate(x);
resvec(end+1) = all(layer.f(x) == y);

% check zonotope
X = zonotope(x,0.01 * eye(4));
Y = layer.evaluate(X);

resvec(end+1) = contains(Y,y);

% gather results
res = all(resvec);


% ------------------------------ END OF CODE ------------------------------
