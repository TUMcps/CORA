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

close all;
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])

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
plot(S,1,'Color',colorblind('r'));
ylabel("plot(S, 1)")

% conHyperplane
i = i + 1;
c = [1, 1];
d = 2;
A = [1, 0];
b = 2.5;
S = conHyperplane(c, d, A, b);

subplot(nrows, ncols, i)
hold on;
xlim([0, 4])
ylim([-1, 1])
plot(S);
title("conHyperplane")

subplot(nrows, ncols, ncols+i)
hold on;
xlim([0, 4])
plot(S,1,'Color',colorblind('r'));

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
plot(S,1,'Color',colorblind('r'));

% conZonotope
i = i + 1;
S = conZonotope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("conZonotope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',colorblind('r'));

% ellipsoid
i = i + 1;
S = ellipsoid(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("ellipsoid")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',colorblind('r'));

% halfspace
i = i + 1;
c = [1, 1];
d = 2;
S = halfspace(c, d);

subplot(nrows, ncols, i)
hold on;
xlim([0, 4])
ylim([-1, 1])
title("halfspace")
plot(S);

subplot(nrows, ncols, ncols+i)
hold on;
xlim([0, 4])
plot(S,1,'Color',colorblind('r'));

% interval
i = i + 1;
S = I;

subplot(nrows, ncols, i)
hold on;
plot(S);
title("interval")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',colorblind('r'));

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
plot(ls,1,'Color',colorblind('r'));

% polytope
i = i + 1;
S = polytope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("polytope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',colorblind('r'));

% polyZonotope
i = i + 1;
S = polyZonotope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("polyZonotope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',colorblind('r'));

% probZonotope
i = i + 1;
Z = [1, 1, -2; 0, 1, 1];
G = [0.6, 1.2; 0.6, -1.2];
S = probZonotope(Z, G);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("probZonotope")
pos = [-57.95, -50.30, 0.76];
set(gca, 'CameraPosition', pos)
xlim([-6, 6])
ylim([-6, 6])

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',colorblind('r'));
set(gca, 'CameraPosition', pos)
xlim([-6, 6])
ylim([-6, 6])

% taylm
i = i + 1;
S = taylm(I);

subplot(nrows, ncols, i)
hold on;
plot(interval(S));
title("taylm")

subplot(nrows, ncols, ncols+i)
hold on;
plot(interval(S),1,'color',colorblind('r'));

% zonoBundle
i = i + 1;

subplot(nrows, ncols, i)
hold on;
xlim([0,3]);
ylim([-1.5,1.5]);

Z1 = zonotope([1.5;0.5], 0.5*eye(2));
Z2 = zonotope([1.5;-0.5], 0.5*eye(2));

plot(Z1, [1, 2], ':',  'Color', colorblind('b'));
plot(Z2, [1, 2], '--',  'Color', colorblind('b'));
plot(zonoBundle({Z1, Z2}), [1, 2], 'Color', colorblind('b'));
title("zonoBundle")

subplot(nrows, ncols, ncols+i)
hold on;
xlim([0,3]);

Z1 = zonotope(1, 0.5);
Z2 = zonotope(2, 0.5);

plot(Z1, 1, ':', 'Color', colorblind('r'));
plot(Z2, 1, '--', 'Color', colorblind('r'));
plot(zonoBundle({Z1,Z2}),1, 'Color', colorblind('r'));

% zonotope
i = i + 1;
S = zonotope(I);

subplot(nrows, ncols, i)
hold on;
plot(S);
title("zonotope")

subplot(nrows, ncols, ncols+i)
hold on;
plot(S,1,'Color',colorblind('r'));

% ------------------------------ END OF CODE ------------------------------
