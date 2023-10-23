function res = testLongDuration_zonotope_minkDiff
% testLongDuration_zonotope_minkDiff - unit test function of minus for
%    approximating the Minkowski difference of two zonotopes or a zonotope
%    with a vector according to [1].
%
% Syntax:
%    res = testLongDuration_zonotope_minkDiff
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Example: 
%    Z1 = zonotope([1 2 2; 0 0 2]);
%    Z2 = zonotope([0 0.5 0.5 0.3; 1 0 0.5 0.2]);
%
%    Z3 = Z1 - Z2;
%    Z4 = Z2 + Z3;
%
%    figure; hold on;
%    plot(Z1,[1 2], 'b');
%    plot(Z2,[1 2], 'r');
%    plot(Z3,[1 2], 'g');
%    plot(Z4,[1 2], 'k');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, conZonotope/minkDiff

% Authors:       Matthias Althoff
% Written:       06-May-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% create zonotopes -- fixed cases in 2D (Minkowski difference is exact in 2D)
Z_m = zonotope([1 1 0 1; 1 0 1 1]);
Z_m_degenerate = zonotope([1 1; 1 0]);
Z_s{1} = zonotope([0 0.5 0; 0 -0.2 0.2]); % see Fig. 2a in [1]
Z_s{2} = zonotope([0 0.5 0; 0 -0.5 0.5]); % see Fig. 2b in [1]
Z_s{3} = zonotope([0 2 0; 0 -0.5 0.5]); % see Fig. 2c in [1]
Z_s{4} = zonotope([0 0.5 0; 0 0 0.5]);
Z_s{5} = zonotope([0 1 0; 0 0 0.5]);
Z_s{6} = zonotope([0 2 0; 0 0 0.5]);
% the following Z_s seem to be incorrect for the MPT toolbox
% Z_s{7} = zonotope([-0.5 0.5 0; 1 0 0.5]); % non-zero center; seems incorrect for MPT toolbox
% Z_s{8} = zonotope([-0.5 1 0; 1 0 0.5]); % non-zero center; seems incorrect for MPT toolbox
% Z_s{9} = zonotope([-0.5 2 0; 1 0 0.5]); % non-zero center; seems incorrect for MPT toolbox
Z_s_degenerate = zonotope([1 0.5; 1 0]);


% convert Z_m to polytope
P_m = polytope(Z_m);
P_m_degenerate = polytope(Z_m_degenerate);

% initialize partial results
resPartial = [];

% define small box
smallBox = zonotope([[0;0],1e-8*eye(2)]);
verySmallBox = zonotope([[0;0],1e-12*eye(2)]);

%% loop through all subtrahends and check for exact result in the 2D case
for iSet = 1:length(Z_s)
    
    % result from polytopes
    % IMPORTANT: Minkowski difference of MPT toolbox seems to be incorrect
    % for minkDiff(Z_s{7},Z_s{9})
    P_res = minkDiff(P_m,polytope(Z_s{iSet}));
    
    % set considered types
    typeSet = {'inner', 'outer'};

    % loop over all types
    for iType = 1:length(typeSet)

        % set type
        type = typeSet{iType};

        % compute result
        Z_res = minkDiff(Z_m, Z_s{iSet}, type);

        % check whether Minkowski difference returns the empty set
        if representsa(Z_res,'emptySet')
            % check if polytope solution is empty as well
            resPartial(end+1) = representsa(P_res,'emptySet');
        else
            % enclosure check (Definition of Minkoswki difference)
            resPartial(end+1) = contains(Z_m + smallBox, Z_res + Z_s{iSet});

            % enclosure check (comparison with polytope solution)
            resPartial(end+1) = contains(P_res + smallBox, polytope(Z_res));

            % enclosure check (comparison with polytope solution; other direction)
            resPartial(end+1) = contains(polytope(Z_res) + smallBox, P_res + verySmallBox);

    %         % for debugging:
    %         figure
    %         hold on
    %         plot(Z_m);
    %         plot(Z_res, [1 2], 'k');
    %         plot(P_res, [1 2], 'r');
    %         plot(Z_s{iSet}, [1 2], 'g');
    %         plot(Z_res + Z_s{iSet}, [1 2], 'k');
        end
    end    
end

%% minuend is a degenerate zonotope and the subtrahend is not
% result from polytopes
P_res = minkDiff(P_m_degenerate,polytope(Z_s{1}));

% set considered types
typeSet = {'inner', 'outer', 'outer:coarse', 'inner:conZonotope', 'inner:RaghuramanKoeln', 'inner:RaghuramanKoeln_weighted', 'inner:incremental'};

% loop over all types
for iType = 1:length(typeSet)
    
    % set type
    type = typeSet{iType};
    
    % compute result
    try 
        Z_res = minkDiff(Z_m_degenerate, Z_s{1}, type);
    catch
        % as expected, this is not possible
        Z_res = [];
    end

    % the result should be empty
    if representsa(Z_res,'emptySet')
        % check if polytope solution is empty as well
        resPartial(end+1) = representsa(P_res,'emptySet');
    end
end

%% minuend and subtrahend are degenerate in 2D (result should be exact)
% result from polytopes
P_res = minkDiff(P_m_degenerate,polytope(Z_s_degenerate));

% set considered types
typeSet = {'inner', 'outer', 'outer:coarse'};

% loop over all types
for iType = 1:length(typeSet)
    % set type
    type = typeSet{iType};
    
    % compute result
    Z_res = minkDiff(Z_m_degenerate, Z_s_degenerate, type);

    % enclosure check (Definition of Minkoswki difference)
    resPartial(end+1) = contains(Z_m_degenerate + smallBox, Z_res + Z_s_degenerate);

    % enclosure check (comparison with polytope solution); not applicable, see
    % above
    resPartial(end+1) = contains(P_res + smallBox, polytope(Z_res));

    % enclosure check (comparison with polytope solution; other direction); not applicable, see
    % above
    resPartial(end+1) = contains(polytope(Z_res) + smallBox, P_res);
end

%% create zonotopes -- fixed case in 3D (no halfspace is redundant; check under-approximation)
% define small box
smallBox = zonotope([[0;0;0],1e-6*eye(3)]);
% create minuend
Z_m = zonotope([zeros(3,1),ones(3,1),eye(3)]);
% create subtrahend
%Z_s = 1/5*zonotope([zeros(3,1),[1; -1; 0],eye(3)]);
%Z_s = 0.228*zonotope([zeros(3,1),[1; -1; 0],eye(3)]); % 0.229 is not working
%Z_s = 0.229*zonotope([zeros(3,1),[1; -1; 0],eye(3)]); % 0.229 is not working
Z_s = 1/4*zonotope([zeros(3,1),[-1; 1; 1],eye(3)]);

% set considered types
typeSet = {'inner', 'inner:conZonotope', 'inner:RaghuramanKoeln', 'inner:RaghuramanKoeln_weighted', 'inner:incremental'};

% loop over all types
for iType = 1:length(typeSet)
    
    % set type
    type = typeSet{iType};

    % compute result
    Z_res = minkDiff(Z_m, Z_s, type);
    
    % enclosure check (Definition of Minkoswki difference)
    resPartial(end+1) = contains(Z_m + smallBox, Z_res + Z_s);
end

% figure
% plot(Z_res,[1 2 3])
% figure
% plot(P_res,[1 2 3])

%% create zonotopes -- random cases in 3D (check under-approximation)
for iSet = 1:10
    % create minuend
    Z_m = zonotope.generateRandom('Dimension',3,'NrGenerators',5);
    
    % create subtrahend
    Z_s = enlarge(zonotope.generateRandom('Dimension',3,'NrGenerators',5), 0.2);
    
    % set considered types
    typeSet = {'inner', 'inner:conZonotope', 'inner:RaghuramanKoeln', 'inner:RaghuramanKoeln_weighted', 'inner:incremental'};

    % loop over all types
    for iType = 1:length(typeSet)
        
        % set type
        type = typeSet{iType};
    
        % compute result
        Z_res = minkDiff(Z_m, Z_s, type);

        % enclosure check (Definition of Minkoswki difference)
        if ~representsa(Z_res,'emptySet')
            resPartial(end+1) = contains(Z_m + smallBox, Z_res + Z_s);
        else
            % if the result is empty, it has to hold that Z_s is not contained
            % in Z_m when the centers are equal
            % This portion is commented out due to errors in the MPT toolbox
            % resPartial(end+1) = ~contains(Z_m, Z_s + smallBox + center(Z_m) - center(Z_s));
        end

        if resPartial(end) ~= 1
            % MPT toolbox often wrongfully returns that a polytope is empty
            path = pathFailedTests(mfilename());
            save(path,'Z_m','Z_s');
            disp(['error in under-approximation of method ',type]);
        end
    end
end

%% create zonotopes -- random cases in 3D (over-approximation)
for iSet = 1:10
    % create minuend
    Z_m = zonotope.generateRandom('Dimension',3,'NrGenerators',4);
    
    % create subtrahend
    Z_s = enlarge(zonotope.generateRandom('Dimension',3,'NrGenerators',4), 0.2);
    
    % compute MPT result
    try % MPT sometimes crashes
        P_res = minkDiff(polytope(Z_m),polytope(Z_s));
    
        % set considered types
        typeSet = {'outer', 'outer:coarse'};

        % loop over all types
        for iType = 1:length(typeSet)
            
            % set type
            type = typeSet{iType};

            % compute result
            Z_res = minkDiff(Z_m, Z_s, type);

            % enclosure check (over-approximation)
            if ~representsa(Z_res,'emptySet')
                resPartial(end+1) = contains(Z_res + smallBox, P_res);
            else
                % if the result is empty, it has to hold that Z_s is not contained
                % in Z_m when the centers are equal
                % This is not always true since the empty set computation
                % of the MPT Toolbox is buggy
                % This portion is commented out due to errors in the MPT toolbox
                % resPartial(end+1) = ~contains(Z_m, Z_s + smallBox + center(Z_m) - center(Z_s));
            end

            if resPartial(end) ~= 1
                % MPT toolbox often wrongfully returns that a polytope is empty
                path = pathFailedTests(mfilename());
                save(path,'Z_m','Z_s');
                disp(['error in over-approximation of method ',type]);
            end
        end
    catch
        disp('MPT Toolbox failed');
    end
end


%result of all tests
res = all(resPartial);

% ------------------------------ END OF CODE ------------------------------
