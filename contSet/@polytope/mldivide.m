function P_out = mldivide(P1, P2)
% mldivide - overloaded '\' operator for the set difference between two
%    polytopes, i.e., the operation
%    P1 \ P2 = { x | x \in P1, x \not\in P2 }
%
% Syntax:
%    P_out = P1 \ P2
%    P_out = mldivide(P1,P2)
%
% Inputs:
%    P1 - polytope object
%    P2 - polytope object
%
% Outputs:
%    P_out - polytope object
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[4;6;4;6;4]);
%
%    P = P1 \ P2
%
% Reference: MPT-Toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Viktor Kotsev
% Written:       09-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{P1,'att','polytope'},...
                {P2,'att','polytope'}});

% polytopes must be in the same dimension
equalDimCheck(P1,P2);

% dimension
n = dim(P1);

% empty set cases
if representsa_(P1,'emptySet')
    % P1 empty -> P1 \ P2 empty
    P_out = polytope.empty(n);
    return
elseif representsa_(P2,'emptySet')
    % P2 empty -> P1 \ P2 = P1
    P_out = P1;
    return
end

% compute minimal H-representation
P1 = compact_(P1,'all',1e-9);
P2 = compact_(P1,'all',1e-9);

% read out constraints
P1_A = P1.A; P1_b = P1.b;
P2_A = P2.A; P2_b = P2.b;
P_out = polytope(zeros(0,n),[]);

% both are fully dimensional, just flip half-spaces and check
% full dimensionality of the results
if isFullDim(P1) && isFullDim(P2)
	for i = 1:length(P2_b)
		A = [-P2_A(i, :); P2_A(1:i-1, :); P1_A];
		b = [-P2_b(i); P2_b(1:i-1); P1_b];
        Pi = polytope(A, b);
		if isFullDim(Pi)
			P_out = Pi;
		end
	end

elseif isFullDim(P1) || (~isempty(P2.Ae) && ~isempty(P1.Ae) && ...
    rank(null(P2.Ae)) < rank(null(P1.Ae)))
    % 1) P1 is full-dim, P2 is low-dim
    % 2) both are low-dim, but P2 has a smaller dimension of
    %    the affine hull
    %
    % unless we support open half-space, the set
    % difference is equal to P1
    P_out = P1;

elseif isFullDim(P2) || (~isempty(P2.Ae) && ~isempty(P1.Ae) && ...
    rank(null(P2.Ae)) > rank(null(P1.Ae)) )
    % 1) P1 is lower-dimensional, but P2 is full dimensional
    % 2) both are lower dimensional, but P2 has higher affine dimension
    %
    % open half-spaces of P2 need to be shifted and emptienies has to be
    % checked 
    P1_Ae = P1.Ae;
    P1_be = P1.be;
    shift_tol = 1e-12;
    for i = 1:length(P1_b)
        A = [-P2_A(i, :); P2_A(1:i-1, :); P1_A];
        b = [-P2_b(i)-shift_tol; P2_b(1:i-1); P1_b];
        if ~isempty([A b]) || ~isempty([P1.Ae P1.be])
            % shift the half-space back
            b(1) = b(1)+shift_tol;
            Pi = polytope(A, b);
            %Pi = polytope('A', A, 'b', b, 'Ae', P1_Ae 'be', P1_be );
            P_out = [P_out Pi];

        end
    end
elseif rank([P2.Ae P2.be; P1.Ae P1.be]) > max(size(P2.Ae, 1), size(P1.Ae, 1))
	% both are lower dimensional, but their affine hulls do not intersect,
	% hence P1 \ P2 = P1
	P_out = P1;

elseif isempty(P2.Ae) && isempty(P1.Ae) && P1<=P2
	% both are lower-dimensional, but with no affine subspace, most
	% probably a vertex
	P_out = polytope();
    
else
	% P1 and P2 are both lower-dimensional and affdim(P1)<=affdim(P2)
	
	% project P1 and P2 on the fully-dimensional null space of P2
	F = null(P2.Ae);
    x0 = P2.Ae\P2.be;
    P1_Z = polytope(P1_A*F, P1_b-P1_A*x0);
	%P1_Z = polytope('A', P1_A*F, 'b' P1_b-P1_A*x0, 'Ae', P1.Ae*F, 'be' P1.be-P1.Ae*x0);
	P2_Z = polytope(P2_A*F, P2_b-P2_A*x0);

	% comptute set difference where P2_Z is now fully dimensional
	P_out = mldivide(P1_Z, P2_Z);
	
	% project back on the original space
	res_new = [];
	for i = 1:numel(P_out)
		if ~isempty(P_out(i))
			A = P_out(i).A;
			b = P_out(i).b;
            Pi = polytope(A*pinv(F), b+A*pinv(F)*x0)
			%Pi = polytope('A', A*pinv(F),'b' b+A*pinv(F)*x0, ...
			%	'Ae', P1.Ae, 'be', P1.be);
			res_new = [res_new Pi];
		end
	end
	P_out = res_new;
end

% ------------------------------ END OF CODE ------------------------------
