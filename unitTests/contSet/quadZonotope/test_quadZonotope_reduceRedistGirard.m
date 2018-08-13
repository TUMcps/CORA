function res = test_quadZonotope_reduceRedistGirard()
% test_quadZonotope_reduceRedistGirard - unit test for zonotope reduction
%                                        with the method "redistGirard"
%
% Syntax:  
%    res = test_quadZonotope_reduceRedistGirard()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:           Niklas Kochdumper
% Written:          17-January-2018
% Last update:      ---
% Last revision:    ---

%------------- BEGIN CODE --------------


    res = true;

    % ANALYTIC TESTS

    options.reduceRedistGirard = 'none';

    % Test 1
    G = [2 0;0 3];
    Grest = [1 -3;2 -1];

    qZ = quadZonotope(zeros(2,1),G,zeros(2,2),zeros(2,1),Grest);
    qZred = reduce(qZ,'redistGirard',0,options);
    c = center(qZred);
    [G_,Gquad_,Gsquare_,Grest_] = generators(qZred);

    G_sol = [5 0;0 6];
    Grest_sol = [1;0];

    if any(c) || any(any(Gquad_)) || any(any(Gsquare_)) || ...
       any(any(G_sol-G_)) || any(any(Grest_sol-Grest_))
        res = false;
        disp('test_quadZonotope_reduceGirard failed');
        return;
    end

    % Test 2
    G = [0 1;2 1];
    Grest = [1 3;-4 0];

    qZ = quadZonotope(zeros(2,1),G,zeros(2,2),zeros(2,1),Grest);
    qZred = reduce(qZ,'redistGirard',0,options);
    c = center(qZred);
    [G_,Gquad_,Gsquare_,Grest_] = generators(qZred);

    G_sol = [0 3;6 3];
    Grest_sol = [-2;2];

    if any(c) || any(any(Gquad_)) || any(any(Gsquare_)) || ...
       any(any(abs(G_sol-G_) > 1e-10)) || any(any(abs(Grest_sol-Grest_) > 1e-10))
        res = false;
        disp('test_quadZonotope_reduceGirard failed');
        return;
    end

end