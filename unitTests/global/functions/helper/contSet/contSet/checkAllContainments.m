function checkAllContainments(S, Sdeg, Sempty, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets)
% checkAllContainments - given a set representation S, checks whether the
%    containment check S' < S can be performed for every other set
%    representation S' for which the check is implemented (only for sets in
%    R^2).
%    The set representations 'affine', 'contSet', 'probZonotope', and 'zoo'
%    are *not* tested.
%    Furthermore, this script does *not* check the formal correctness of
%    the 'approx' method, or that of any algorithm in additionalAlgorithms;
%    it only verifies that they don't throw any MATLAB/CORA errors.
%    It also does *not* check the optional 'scaling' output, though it does
%    check the 'cert' output.
%    Additionally, it does *not* check additionalAlgorithms on
%    point-containment.
%    Also, note that taylm has been disabled for now - too many problems.
%
% Syntax:
%    checkAllContainments(S, Sdeg, Sempty, containsNonDeg, ...
%                         implementedSets, setsNonExact, ...
%                         additionalAlgorithms, ...
%                         additionalAlgorithms_specificSets)
%
% Inputs:
%    S - Basic set of the given type; should more or less be as big as a
%        unit hypercube.
%        In the case of non-compact sets: it should contain the negative
%        x_1-axis, but *not* the positive x_2-axis.
%    Sdeg - Same as S, but with a degenerate version of the set (or an
%           empty object if the set can not be degenerate, by
%           construction).
%           In the case of non-compact sets: it should contain the negative
%           x_1-axis, but *not* the positive x_2-axis.
%    Sempty - Empty instance of the set
%    implementedSets - cell array of character arrays describing all set
%                      representations S' for which the containment S' < S
%                      is implemented (beyond point-containment, emptyset,
%                      and fullspace)
%    setsNonExact - cell array of character arrays describing all set
%                   representations from implementedSets for which there is
%                   no 'exact' algorithm (or, at least, there might be some
%                   cases where it fails to be exact).
%                   Note that one can add 'point' here, in case
%                   point-containment is not exact.
%                   For these set representations, the correctness of the
%                   containment is **not** checked; only the correct
%                   execution of the algorithm.
%    additionalAlgorithms - cell array of character arrays describing all
%                           implemented algorithms, beyond 'exact' and
%                           'approx'
%    additionalAlgorithms_specificSets - For each entry in
%                   additionalAlgorithms, you may specify on which set
%                   representations to apply the algorithm in question. If
%                   you do not specify anything, all set representations
%                   will be tested.
%                   See unitTests/zonotope/testLong_zonotope_contains for
%                   an example of how this is used.
%
% Outputs:
%    -
%
% Example:
%    %% For an example use of checkAllContainments, we refer to
%    % unitTests/conZontope/testLong_conZonotope_contains
%    % as it requires the entire code showing off every aspect of
%    % checkAllContainments
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: unitTests/conZonotope/testLong_conZonotope_contains

% Authors:       Adrian Kulmburg
% Written:       12-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Set global tolerance for every containment check
tol = 1e-5;

%% Definitions of all tested sets

% Points
point_contained = [0;0]; % Origin should be contained in basic set
point_notContained = [100;0]; % A far away point that is not contained

% Other compact sets
% List of 'small' sets that are certainly contained within the basic set
C_small = capsule([0;0],[0.01;0],0.01);
I_small = interval([-0.01;-0.01],[0.01;0.01]);
cPZ_small = conPolyZono(I_small);
cZ_small = conZonotope(I_small);
E_small = ellipsoid(0.001*eye(2));
syms x y; eqs = x^2 + y^2 - 0.001; ls_small = levelSet(eqs,[x;y],'<=');
P_small = 0.01*polytope([1 1 -1 -1;1 -1 1 -1]);
pZ_small = polyZonotope(I_small);
SpS_small = spectraShadow(I_small);
%tay_small = taylm(I_small);
zB_small = zonoBundle(I_small);
Z_small = zonotope(I_small);

smallSets = {C_small, I_small, cPZ_small, cZ_small, E_small, ls_small,...
                P_small, pZ_small, SpS_small,...
                ...tay_small,...
                zB_small, Z_small};

% List of 'big' sets that are certainly not contained within the basic set
C_big = capsule([0;0],[100;0],100);
I_big = interval([-100;-100],[100;100]);
cPZ_big = conPolyZono(I_big);
cZ_big = conZonotope(I_big);
E_big = ellipsoid(100*eye(2));
syms x y; eqs = x^2 + y^2 - 100; ls_big = levelSet(eqs,[x;y],'<=');
P_big = 100*polytope([1 1 -1 -1;1 -1 1 -1]);
pZ_big = polyZonotope(I_big);
SpS_big = spectraShadow(I_big);
%tay_big = taylm(I_big);
zB_big = zonoBundle(I_big);
Z_big = zonotope(I_big);

bigSets = {C_big, I_big, cPZ_big, cZ_big, E_big, ls_big,...
                P_big, pZ_big, SpS_big,...
                ...tay_big,...
                zB_big, Z_big};

% List of 'small' **degenerate** sets that are certainly contained
C_degSmall = capsule.empty(2); % capsules can not be degenerate otherwise
I_degSmall = interval([-0.01;0],[0.01;0]);
cPZ_degSmall = conPolyZono(I_degSmall);
cZ_degSmall = conZonotope(I_degSmall);
E_degSmall = [1 0; 0 0] * ellipsoid(0.001*eye(2));
syms x y; eqs = [x^2 + y^2 - 0.001; y; -y]; ls_degSmall = levelSet(eqs,[x;y],'<=');
P_degSmall = 0.01*polytope([-1 1; 0 0]);
pZ_degSmall = polyZonotope(I_degSmall);
SpS_degSmall = spectraShadow(I_degSmall);
%tay_degSmall = taylm(I_degSmall);
zB_degSmall = zonoBundle(I_degSmall);
Z_degSmall = zonotope(I_degSmall);

degSmallSets = {C_degSmall, I_degSmall, cPZ_degSmall, cZ_degSmall, ...
                E_degSmall, ls_degSmall, P_degSmall, pZ_degSmall, ...
                SpS_degSmall,...
                ...tay_degSmall,...
                zB_degSmall, Z_degSmall};

% List of 'big' **degenerate** sets that are certainly *not* contained
% C_degBig = --- ; % capsules can not be degenerate 
I_degBig = interval([-100;0],[100;0]);
cPZ_degBig = conPolyZono(I_degBig);
cZ_degBig = conZonotope(I_degBig);
E_degBig = [1 0; 0 0] * ellipsoid(100*eye(2));
syms x y; eqs = [x^2 + y^2 - 100; y; -y]; ls_degBig = levelSet(eqs,[x;y],'<=');
P_degBig = 100*polytope([-1 1; 0 0]);
pZ_degBig = polyZonotope(I_degBig);
SpS_degBig = spectraShadow(I_degBig);
%tay_degBig = taylm(I_degBig);
zB_degBig = zonoBundle(I_degBig);
Z_degBig = zonotope(I_degBig);

degBigSets = {... %C_degBig,
                I_degBig, cPZ_degBig, cZ_degBig, ...
                E_degBig, ls_degBig, P_degBig, pZ_degBig, ...
                SpS_degBig,...
                ...tay_degBig,...
                zB_degBig, Z_degBig};


% List of empty instances of all sets
C_empty = capsule.empty(2);
I_empty = interval.empty(2);
cPZ_empty = conPolyZono.empty(2);
cZ_empty = conZonotope.empty(2);
E_empty = ellipsoid.empty(2);
ls_empty = levelSet.empty(2);
P_empty = polytope.empty(2);
pZ_empty = polyZonotope.empty(2);
SpS_empty = spectraShadow.empty(2);
%tay_empty = taylm.empty(2);
zB_empty = zonoBundle.empty(2);
Z_empty = zonotope.empty(2);

emptySets = {C_empty, I_empty, cPZ_empty, cZ_empty, E_empty, ls_empty,...
                P_empty, pZ_empty, SpS_empty,...
                ...tay_empty,...
                zB_empty, Z_empty};

% We can now perform all the actual containment checks

%% Empty set and fullspace containment
% Empty set containment checks
aux_executeContainmentChecks_emptySet(S, Sdeg, Sempty, additionalAlgorithms);

% Fullspace containment checks
aux_executeContainmentChecks_fullspace(S, Sdeg, Sempty, additionalAlgorithms)

%% Point-containment
% Check for exact containment (if implemented)
exactIsImplemented = false;
try
    % Origin must be contained, unless S is empty
    res = contains(S,point_contained,'exact',tol);
    [res_cert, cert] = contains(S,point_contained,'exact',tol);
    if representsa(S, 'emptySet')
        assert(~res);
        assert(~res_cert);
        assert(cert);
    else
        assert(res);
        assert(res_cert);
        assert(cert);
    end

    % Far away point must not be contained, except if S is the fullspace of
    % course
    res = contains(S,point_notContained,'exact',tol);
    [res_cert, cert] = contains(S,point_notContained,'exact',tol);
    if representsa(S,'fullspace')
        assert(res);
        assert(res_cert);
    else
        assert(~res);
        assert(~res_cert);
    end
    assert(cert);

    % If we came that far, it means that an exact algorithm for
    % point-containment must have been implemented; check that this is
    % intended:
    exactIsImplemented = true;
catch ME
    % The code above is only allowed to fail if there is no exact algorithm
    % for point-containment
    if ~any(strcmp(setsNonExact, 'point'))
        rethrow(ME);
    end
end

% If the point-containment behaves like it is exact, make sure that this is
% intended:
% @@@
% DO NOT REMOVE THIS CHECK!!!
% If your contains-method fails this check, then *your* method is faulty,
% not the present check; if your method can't recognize that a
% super-far-away point is not contained, then it is worthless.
% @@@
assert(~(exactIsImplemented && any(strcmp(setsNonExact, 'point'))))

% There **MUST** be at least an approximative way of checking for
% point-containment. Otherwise, the contains function is useless.
% Try it out, without necessarily checking for the correctness of the
% result:
contains(S,point_contained,'approx');
[~,~] = contains(S,point_contained,'approx');

%% Basic containment checks
% If S is non-degenerate, it must contain smallSets, otherwise not
% (Currently, the only case where S could be degenerate here is when S is a
% conHyperplane)
shouldBeContained = isFullDim(S);
aux_executeContainmentChecks_basic(S, smallSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
% S must contain degSmallSets, emptySets
shouldBeContained = true;
aux_executeContainmentChecks_basic(S, degSmallSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
aux_executeContainmentChecks_basic(S, emptySets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
% S must *not* contain bigSets, degBigSets
shouldBeContained = false;
aux_executeContainmentChecks_basic(S, bigSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
aux_executeContainmentChecks_basic(S, degBigSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

% Sdeg must contain emptySets
shouldBeContained = true;
aux_executeContainmentChecks_basic(Sdeg, emptySets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
% If Sdeg is not empty, it must contain smallDegSets, otherwise not
shouldBeContained = ~representsa(Sdeg, 'emptySet');
aux_executeContainmentChecks_basic(Sdeg, degSmallSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
% Sdeg must *not* contain smallSets, bigSets, degBigSets
shouldBeContained = false;
aux_executeContainmentChecks_basic(Sdeg, smallSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
aux_executeContainmentChecks_basic(Sdeg, bigSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
aux_executeContainmentChecks_basic(Sdeg, degBigSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

% Sempty must contain emptySets
shouldBeContained = true;
aux_executeContainmentChecks_basic(Sempty, emptySets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
% Sempty must *not* contain smallSets, smallDegSets, bigSets, degBigSets,
% unless they are empty (this is checked internally)
shouldBeContained = false;
aux_executeContainmentChecks_basic(Sempty, smallSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
aux_executeContainmentChecks_basic(Sempty, degSmallSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
aux_executeContainmentChecks_basic(Sempty, bigSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);
aux_executeContainmentChecks_basic(Sempty, degBigSets, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets);

end


% Auxiliary functions -----------------------------------------------------

function aux_executeContainmentChecks_emptySet(S, Sdeg, Sempty, additionalAlgorithms)
    % Checks if S, Sdeg, Sempty contain the empty set, and wether S and
    % Sdeg are *not* contained in the empty set, but Sempty is contained in
    % the empty set.
    % Does this for 'exact', 'approx', and every algorithm in
    % additionalAlgorithms.

    es = emptySet(2);

    algorithms = ['exact' 'approx' additionalAlgorithms];

    for i=1:length(algorithms)
        algo = algorithms{i};

        % No matter the algorithm, the empty set needs to be contained
        S_res = contains(S,es,algo);
        [S_res_cert, S_cert] = contains(S,es,algo);
        Sdeg_res = contains(Sdeg,es,algo);
        [Sdeg_res_cert, Sdeg_cert] = contains(S,es,algo);
        Sempty_res = contains(Sempty,es,algo);
        [Sempty_res_cert, Sempty_cert] = contains(Sempty,es,algo);

        res = S_res && S_res_cert && Sdeg_res && Sdeg_res_cert && Sempty_res && Sempty_res_cert;
        cert = S_cert && Sdeg_cert && Sempty_cert;

        assert(res);
        assert(cert);
    end

    % Now, we check that S is not contained in the empty set, but
    % that Sempty is contained; for Sdeg, it depends whether it is empty or
    % not.
    % We also check the scaling; might as well
    % Starting with the exact algorithm
    S_res = contains(es,S);
    [S_res_cert, S_cert] = contains(es,S);
    [S_res_scaling, S_cert_scaling, S_scaling] = contains(es,S);
    Sdeg_res = contains(es,Sdeg);
    [Sdeg_res_cert, Sdeg_cert] = contains(es,Sdeg);
    [Sdeg_res_scaling, Sdeg_cert_scaling, Sdeg_scaling] = contains(es,Sdeg);
    Sempty_res = contains(es,Sempty);
    [Sempty_res_cert, Sempty_cert] = contains(es,Sempty);
    [Sempty_res_scaling, Sempty_cert_scaling, Sempty_scaling] = contains(es,Sempty);

    degIsEmpty = representsa(Sdeg, 'emptySet');
    % If degIsempty is true, Sdeg may be contained in the empty set

    % If S is empty itself (might only happen if S = emptySet), we can skip
    % the first check
    if representsa(S, 'emptySet')
        res = true;
    else
        res = ~S_res && ~S_res_cert && ~S_res_scaling;
    end
    res = res && (Sdeg_res==degIsEmpty) && (Sdeg_res_cert==degIsEmpty) && (Sdeg_res_scaling==degIsEmpty);
    res = res && Sempty_res && Sempty_res_cert && Sempty_res_scaling;

    cert = S_cert && S_cert_scaling && Sdeg_cert && Sdeg_cert_scaling && Sempty_cert && Sempty_cert_scaling;

    % If S is the empty set, we need a special check
    if representsa(S, 'emptySet')
        scaling = (S_scaling == 0);
    else
        scaling = (S_scaling == inf); % S cannot be contained in any scaling of the empty set
    end
    if degIsEmpty
        scaling = scaling && (Sdeg_scaling == 0); % If Sdeg is empty, scaling should be 0
    else
        scaling = scaling && (Sdeg_scaling == inf); % Otherwise, same as S
    end
    scaling = scaling && (Sempty_scaling == 0);

    assert(res);
    assert(cert);
    assert(scaling);
end

function aux_executeContainmentChecks_fullspace(S, Sdeg, Sempty, additionalAlgorithms)
    % Checks if S, Sdeg, Sempty do *not*  contain the full space, and
    % wether S, Sdeg, Sempty are contained in the fullspace.
    % Does this for 'exact', 'approx', and every algorithm in
    % additionalAlgorithms.

    tol = 1e-6;

    fs = fullspace(2);

    algorithms = ['exact' 'approx' additionalAlgorithms];

    for i=1:length(algorithms)
        algo = algorithms{i};

        % No matter the algorithm, the full space can *not* be contained
        % (except, of course, if it is the fullspace)
        S_res = contains(S,fs,algo);
        [S_res_cert, S_cert] = contains(S,fs,algo);
        Sdeg_res = contains(Sdeg,fs,algo);
        [Sdeg_res_cert, Sdeg_cert] = contains(Sdeg,fs,algo);
        Sempty_res = contains(Sempty,fs,algo);
        [Sempty_res_cert, Sempty_cert] = contains(Sempty,fs,algo);

        if representsa(S, 'fullspace')
            res = S_res && S_res_cert && ~Sdeg_res && ~Sdeg_res_cert && ~Sempty_res && ~Sempty_res_cert;
        else
            res = ~S_res && ~S_res_cert && ~Sdeg_res && ~Sdeg_res_cert && ~Sempty_res && ~Sempty_res_cert;
        end
        cert = S_cert && Sdeg_cert && Sempty_cert;

        assert(res);
        assert(cert);
    end

    % Now, we check that S, Sdeg, Sempty are contained in the full space.
    % We also check the scaling; might as well
    % Starting with the exact algorithm
    S_res = contains(fs,S);
    [S_res_cert, S_cert] = contains(fs,S);
    [S_res_scaling, S_cert_scaling, S_scaling] = contains(fs,S);
    Sdeg_res = contains(fs,Sdeg);
    [Sdeg_res_cert, Sdeg_cert] = contains(fs,Sdeg);
    [Sdeg_res_scaling, Sdeg_cert_scaling, Sdeg_scaling] = contains(fs,Sdeg);
    Sempty_res = contains(fs,Sempty);
    [Sempty_res_cert, Sempty_cert] = contains(fs,Sempty);
    [Sempty_res_scaling, Sempty_cert_scaling, Sempty_scaling] = contains(fs,Sempty);

    res = S_res && S_res_cert && S_res_scaling;
    res = res && Sdeg_res && Sdeg_res_cert && Sdeg_res_scaling;
    res = res && Sempty_res && Sempty_res_cert && Sempty_res_scaling;

    cert = S_cert && S_cert_scaling && Sdeg_cert && Sdeg_cert_scaling && Sempty_cert && Sempty_cert_scaling;

    scaling = (S_scaling == 0); % S is contained in any scaling of fullspace,
                                % so the infimum over all these values is
                                % 0. A similar argument works for Sdeg and
                                % Sempty
    scaling = scaling && (Sdeg_scaling == 0); 
    scaling = scaling && (Sempty_scaling == 0);

    assert(res);
    assert(cert);
    assert(scaling);
end

function aux_executeContainmentChecks_basic(basicSet, setCollection, shouldBeContained, implementedSets, setsNonExact, additionalAlgorithms, additionalAlgorithms_specificSets)
    % Checks if basicSet contains (or not, depending on shouldBeContained)
    % all sets from setCollection, if they are implemented

    tol = 1e-6;
    N = 3; % Very low number of iterations for opt and sampling methods

    for i = 1:length(setCollection)
        set = setCollection{i};
    
        exactIsImplemented = false;
        try
            res = contains(basicSet,set,'exact',tol);
            [res_cert, cert] = contains(basicSet,set,'exact',tol);

            if representsa(set,'emptySet')
                % Empty sets must always be contained
                assert(res == true);
                assert(res_cert == true);
                assert(cert);
            elseif representsa(basicSet,'emptySet')
                % The empty set contains nothing, except the empty set,
                % which is treated above anyhow
                assert(res == false);
                assert(res_cert == false);
                assert(cert);
            elseif representsa(basicSet,'fullspace')
                % The fullspace contains everything
                assert(res == true);
                assert(res_cert == true);
                assert(cert);
            else
                assert(res == shouldBeContained);
                assert(res_cert == shouldBeContained);
                assert(cert);
            end
    
            % Also, check if the 'approx' method works
            contains(basicSet,set,'approx');
            [~,~] = contains(basicSet,set,'approx');
    
            % If we came that far, it means that an exact algorithm for
            % containment must have been implemented; check that this is
            % intended:
            exactIsImplemented = true;
        
        catch ME
            % The containment check may only fail if:
            if ~any(strcmp(implementedSets, class(set)))
                % 1) the containment check is not implemented at all -> ignore
                continue
            elseif any(strcmp(setsNonExact, class(set)))
                % 2) the exact containment check is not available -> make sure
                % that at least the 'approx' algorithm works
                contains(basicSet,set,'approx');
                [~,~] = contains(basicSet,set,'approx');
            else
                rethrow(ME);
            end
        end

        % Also need to verify that the set is non-empty - if it is, the
        % method may have a direct way of checking for containment in this
        % case. Same holds for cases where the circumbody is empty
        isEmpty = representsa(set,'emptySet') | representsa(basicSet, 'emptySet');
    
        % If the containment behaves like it is exact, make sure that this is
        % intended:
        % @@@
        % DO NOT REMOVE THIS CHECK!!!
        % If your contains-method fails this check, then *your* method is
        % faulty, not the present check; if your method can't recognize that a
        % super-far-away set is not contained, then it is worthless.
        % @@@
        if ~isEmpty
            assert(~(exactIsImplemented && any(strcmp(setsNonExact, class(set)))))
        else
            % Empty sets must always be recognized
            assert(exactIsImplemented);
        end
    end

    % It only remains to check the additional algorithms; note we are not
    % interested here whether they deliver the correct result, we only
    % check correct execution
    for i=1:length(additionalAlgorithms)
        algorithm = additionalAlgorithms{i};
        algorithm_implementedSets = additionalAlgorithms_specificSets{i};

        for j = 1:length(setCollection)
            set = setCollection{j};

            algorithmIsImplemented = false;
            try
                res = contains(basicSet,set,algorithm,tol,N);
                [res_cert, cert] = contains(basicSet,set,algorithm,tol,N);

                % If we came that far, it means that an exact algorithm for
                % point-containment must have been implemented; check that this is
                % intended:
                algorithmIsImplemented = true;

            catch ME
                % The containment check may only fail if it is not
                % implemented
                if ~any(strcmp(algorithm_implementedSets, class(set))) && ~isempty(algorithm_implementedSets)
                    continue
                else
                    rethrow(ME);
                end
            end

            % Also need to verify that the set is non-empty - if it is, the
            % method may have a direct way of checking for containment in this
            % case. Same holds for cases where the circumbody is empty
            isEmpty = representsa(set,'emptySet') | representsa(basicSet, 'emptySet');

            % If the containment behaves like it is implemented, make sure that this is
            % intended:
            % @@@
            % DO NOT REMOVE THIS CHECK!!!
            % If your contains-method fails this check, then *your* method is
            % faulty, not the present check; it means that you indicated
            % the method should not be implemented for the current set
            % representation, yet it does function in practice
            % @@@
            checkToFail = algorithmIsImplemented && any(strcmp(algorithm_implementedSets, class(set))) && isempty(algorithm_implementedSets);
            if ~isEmpty
                assert(~checkToFail);
            else
                % Empty sets must always be recognized
                assert(algorithmIsImplemented);
            end
        end

    end
end


% ------------------------------ END OF CODE ------------------------------
