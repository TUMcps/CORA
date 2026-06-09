function res = example_zonotope_reduceUnderApprox
% example_zonotope_reduceUnderApprox - example of reduceUnderApprox
%
% Syntax:
%    res = example_zonotope_reduceUnderApprox
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
% See also: none

% Authors:       Laura Luetzow
% Written:       05-March-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% methods that shall be evaluated
methodsRed = [...
    "yang","raghuraman","kochdumper","sadraddini", ...
    "cluster","boxCone","scale","nlp",...
   ]; 

% settings for plots
plot_colors = lines;
line_colors = dictionary(["yang","raghuraman","kochdumper","sadraddini", ...
    "cluster","boxCone","scale","nlp"], ...
    {plot_colors(1,:), plot_colors(2,:), plot_colors(3,:), plot_colors(4,:), ...
    plot_colors(1,:), plot_colors(2,:), plot_colors(3,:), plot_colors(4,:)});
line_types = dictionary(["yang","raghuraman","kochdumper","sadraddini", ...
    "cluster","boxCone","scale","nlp"], ...
    {':', ':', ':', ':', '-', '-', '-', '-'});

% create pairs [dim, order] that shall be evaluated
dim = 2;
order = 5;

% order of the reduced zonotope (if negative: add to initial order of
% random zonotope)
orderRed0 = 1; % or 2


% maximum dimension for which the exact zonotope volume is computed
maxDim = 5; %TO-DO: increase!

% maximum generator number for which the exact zonotope volume is computed
maxGenerators = 150;

% distributions for sampling the random zonotope ('uniform', 'exp',
% 'gamma') -> all specified distributions are used for about the same
% amount of zonotopes
dist = ["uniform"]; % or  "exp", "gamma"

fprintf('#####################################################\n');
fprintf('#### Parameters: orderReduce=%d\n', orderRed0);

if orderRed0 < 0 % if negative: add to initial order of random zonotope
    orderRed = order + orderRed0;
else
    orderRed = orderRed0;
end

fprintf('## Parameters: dim=%d, order=%d\n', dim, order);

T = NaN*ones(length(methodsRed),1);
R = NaN*ones(length(methodsRed),1);
Zreds = cell(length(methodsRed),1);

%generate random zonotope using the distribution dist
Zrand = zonotope.generateRandom('Dimension',dim,'NrGenerators',order*dim, "distribution", dist);

%compute volume of zonotope
V_Zrand = aux_volume(Zrand,maxDim,maxGenerators);

for i_method = 1:length(methodsRed)
    method = methodsRed(i_method);

    % compute underapparoximation
    Timer1 = tic;
    Zred = reduceUnderApprox(Zrand, method, orderRed);
    T(i_method) = toc(Timer1);

    % compute volume
    V_red = aux_volume(Zred,maxDim,maxGenerators);
    R(i_method) = aux_VolRatio(V_red,V_Zrand,dim);
    Zreds{i_method} = Zred;
end

% plot specific zonotopes and their underapproximations
aux_plotZonotopes(Zreds, Zrand, line_colors, line_types, order, orderRed0, methodsRed)

% test completed
res = true;
end


% Auxiliary functions -----------------------------------------------------

function [ratio]=aux_VolRatio(Vred,Vrand,dim)
% compute volume ratio

ratio=(Vred/Vrand)^(1/dim);
end


function V=aux_volume(Z,maxDim,maxGenerators)
% compute zonotope volume

if (dim(Z) > maxDim || size(generators(Z),2) > maxGenerators) % && (size(generators(Z),2) / dim(Z) > 1)
    % use interval norm
    V = 2^dim(Z)*prod(sum(abs(generators(Z)),2));
else
    V=volume_(Z, 'exact');
end
end


function aux_plotZonotopes(Zreds, Zrand, line_colors, line_types, order, orderRed, methodsRed)
% plot specific zonotopes and their underapproximations

figure;
hold on;
plot(Zrand,[1 2],'Color', [0.7 0.7 0.7], 'DisplayName', 'Random Zonotope', 'LineWidth', 4);

for i = 1:length(Zreds)
    method = methodsRed(i);
    if sum(abs(Zreds{i}.generators),'all') < 1e-6
        c=Zreds{i}.center;
        plot(c(1),c(2),'*','Color', line_colors{method}, 'DisplayName', method, 'MarkerSize', 6);
    else
        if line_types{method} == ':'
            lw = 1.5;
        else
            lw = 2;
        end
        plot(Zreds{i},[1 2], line_types{method},'Color', line_colors{method}, 'DisplayName', method, 'LineWidth', lw);
    end
end
legend
title(sprintf("Zonotopes Reduced from Order %d to Order %d", order, orderRed))

end


% ------------------------------ END OF CODE ------------------------------
