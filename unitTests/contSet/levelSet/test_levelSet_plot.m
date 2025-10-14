function res = test_levelSet_plot
% test_levelSet_plot - unit test function of plot
%
% Syntax:
%    res = test_levelSet_plot
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
% See also: plot.m

% Authors:       Niklas Kochdumper
% Written:       28-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    tol = 0.01;

    % 2D Equality + Inequality (fixed) ------------------------------------
    
    % construct level set
    syms x1 x2;
    
    eq = [x1^2 + x2^2 - 1; ...
          x2 - 0.5];
    
    ls = levelSet(eq,[x1;x2],{'==';'<='});

    % plot level set
    figure; hold on; box on;
    xlim([-1,1]); ylim([-1,1]);
    plot(ls,[1,2]);

    % extract data from plot
    h = findobj(gca,'Type','line');
    
    % check for correctness
    for i = 1:length(h)
        p = [h(i).XData;h(i).YData];
        assert(all(contains(ls,p,'exact',tol)));
    end

    close


    % 3D Equality Constraints (fixed) -------------------------------------

    % define level sets to test
    syms x1 x2 x3;

    eqs = {};
    eqs{end+1} = [x1^2 + x2^2 - (x3+2)^2; x1^2 + x2^2 - 1.2];
    eqs{end+1} = x1^2 + x2^2 + (x3+1.8)^2 - 1;
    eqs{end+1} = x1^3*x2^3 - (1+0.1*x3^3);

    % loop over all level sets
    for i = 1:length(eqs)

        % construct level set
        if length(eqs{i}) == 1
            ls = levelSet(eqs{i},[x1;x2;x3],'==');
        else
            ls = levelSet(eqs{i},[x1;x2;x3],{'==','<='});
        end

        % plot level set
        figure; hold on; box on;
        xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
        plot(ls,[1,2,3]);

        % extract data from plot
        h = findobj(gca,'Type','patch'); 

        % check correctness
        for j = 1:length(h)
            assert(all(contains(ls,h(j).Vertices','exact',20*tol)))
        end

        close;
    end


    % 3D Equality Constraints (multiple) ----------------------------------

    % construct level set
    syms x1 x2 x3;

    %eq1 = x1^3*x2^3 - (1+0.1*x3^3);   (takes quite long)
    %eq2 = x3 + 0.9;                   (takes quiet long)

    eq1 = x1^2 + x2^2 - (x3+2.2)^2;
    eq2 = x1^3*x2^3 - (0.9+0.1*x3^3);

    ls = levelSet(eq1,[x1;x2;x3],'==') & levelSet(eq2,[x1;x2;x3],'==');

    % plot level set
    figure; hold on; box on;
    xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
    plot(ls,[1,2,3],'r','LineWidth',2);

    % extract data from plot
    h = findobj(gca,'Type','line');
    
    % check for correctness
    for i = 1:length(h)
        p = [h(i).XData;h(i).YData;h(i).ZData];
        assert(all(contains(ls,p,'exact',tol)));
    end

    close


    % 2D Inequality Constraints (fixed) -----------------------------------
    
    % construct level set
    syms x1 x2;
    
    eq = 0.873*(x1 + 0.295)*(x2 - 0.0606)^2 + 0.0034*(x1 + 0.295)^2* ...
         (x2 - 0.0606) - 0.131*(x1 + 0.295)*(x2 - 0.0606)^3 - ...
         0.108*(x1 + 0.295)^3*(x2 - 0.0606) - 0.395*(x1 + 0.295)^3 - ...
         0.727*(x1 + 0.295)^4 + 0.789*(x1 + 0.295)^2*(x2 - 0.0606)^2 + ...
         0.104*(x2 - 0.0606)^3 + 0.256*(x2 - 0.0606)^4 + ...
         0.34*(x1 + 0.295)*(x2 - 0.0606) - 0.00392;
    
    ls = levelSet(eq,[x1;x2],'<=');
    
    % plot level set
    figure; hold on;
    xlim([-1.5,1.5]); ylim([-1.5,1.5]);
    plot(ls);
    
    % extract filled area from the plot
    h = findobj(gca,'Type','patch'); 
    pgon = [];
    
    for i = 1:length(h)
        pgon = pgon | polygon(h(i).Vertices(:,1),h(i).Vertices(:,2));
    end
    
    % check for correctness
    pIn = randPoint(pgon,100);
    pOut = randPoint(subtract(polygon(interval(pgon)),pgon),100);
    
    assert(all(contains(ls,pIn,'exact',tol)));
    assert(all(contains(not(ls),pOut,'exact',tol)));

    close;

    % test completed
    res = true;

% ------------------------------ END OF CODE ------------------------------
