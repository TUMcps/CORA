function example_plot1Dsets()
% example_plot1Dsets - shows all plots possible with 1D sets
%
% Syntax:
%    completed = example_plot1Dsets()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       31-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

figure;

nrows = 2;
ncols = 14;
i = 0;

I = interval([1; 0], [2; 0]);

% capsule
i = i + 1;
S = capsule(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("capsule")
ylabel("plot(S_{1D}, [1, 2])")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));
ylabel("plot(S, 1)")

% conPolyZono
i = i + 1;
S = conPolyZono(I);

subplot(nrows, ncols, i)
hold on;
ylim([-1, 1])
plot(S);
title("conPolyZono")

subplot(nrows, ncols, ncols+i)
hold on;
ylim([-1, 1])
plot(S,1,'Color',CORAcolor('CORA:red'));

% conZonotope
i = i + 1;
S = conZonotope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("conZonotope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));

% ellipsoid
i = i + 1;
S = ellipsoid(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("ellipsoid")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));

% interval
i = i + 1;
S = I;

subplot(nrows, ncols, i)
hold on;
plot(S);
title("interval")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));

% levelset
i = i + 1;

syms x y
eq = (x-1.5)^2 - 0.25; % [1;2]
ls = levelSet(eq,[x;y],"==");

subplot(nrows, ncols, i); 
xlim([0,3]); ylim([-1,1]);

plot(ls);
title("levelset")

subplot(nrows, ncols, ncols+i)
xlim([0,3]); ylim([-1,1]);
plot(ls,1,'Color',CORAcolor('CORA:red'));

% polytope
i = i + 1;
S = polytope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("polytope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));

% polyZonotope
i = i + 1;
S = polyZonotope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("polyZonotope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));

% probZonotope ---
i = i + 1;
Z = [1, 1, -2; 0, 1, 1];
G = [0.6, 1.2; 0.6, -1.2];
S = probZonotope(Z, G);

% plot 1d
subplot(nrows, ncols, i); hold on;
plot(S);
title("probZonotope")
pos = [-57.95, -50.30, 0.76];
set(gca, 'CameraPosition', pos)
xlim([-6, 6]); ylim([-6, 6]);

% plot 2d
subplot(nrows, ncols, ncols+i); hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));
set(gca, 'CameraPosition', pos)
xlim([-6, 6]); ylim([-6, 6])

% taylm
i = i + 1;
S = taylm(I);

subplot(nrows, ncols, i)
hold on;
plot(interval(S));
title("taylm")

subplot(nrows, ncols, ncols+i)
hold on;
plot(interval(S),1,'Color',CORAcolor('CORA:red'));

% zonoBundle
i = i + 1;

% 2d plot
subplot(nrows, ncols, i)
hold on;
xlim([0,3]);
ylim([-1.5,1.5]);

Z1 = zonotope([1.5;0.5], 0.5*eye(2));
Z2 = zonotope([1.5;-0.5], 0.5*eye(2));

plot(Z1, [1, 2], ':',  'Color', CORAcolor('CORA:blue'));
plot(Z2, [1, 2], '--',  'Color', CORAcolor('CORA:blue'));
plot(zonoBundle({Z1, Z2}), [1, 2], 'Color', CORAcolor('CORA:blue'));
title("zonoBundle")

% 1d plot
subplot(nrows, ncols, ncols+i)
hold on;
xlim([0,3]);

Z1 = zonotope(1, 0.5);
Z2 = zonotope(2, 0.5);

plot(Z1, 1, ':', 'Color', CORAcolor('CORA:red'));
plot(Z2, 1, '--', 'Color', CORAcolor('CORA:red'));
plot(zonoBundle({Z1,Z2}),1, 'Color', CORAcolor('CORA:red'));

% zonotope
i = i + 1;
S = zonotope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("zonotope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',CORAcolor('CORA:red'));

% ------------------------------ END OF CODE ------------------------------
