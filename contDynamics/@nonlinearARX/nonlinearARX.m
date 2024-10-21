classdef nonlinearARX < contDynamics
% nonlinearARX class (time-discrete nonlinear ARX model)
%
% Syntax:
%    % only dynamic equation
%    nlnARX = nonlinearARX(fun,dt,dim_y,dim_u,n_p)
%    nlnARX = nonlinearARX(name,fun,dt,dim_y,dim_u,n_p)
%
% Description:
%    Generates a discrete-time nonlinear ARX object (NARX) according  
%    to the following equation:
%       y(k) = f(y(k-1),...,y(k-n_y),u(k),...,u(k-n_u),e(k-1),...,e(k-n_e)) 
%               + e(k)
%
% Inputs:
%    fun    - function handle for the NARX equation with arguments (y,u)
%               y = [y(k-1); ...; y(k-n_p)]: array dim_y x n_p
%               u = [u(k); ...; u(k-n_p)]: array dim_u x (n_p+1)        
%    name   - name of the model
%    dt     - sampling time
%    dim_y  - dimension of the output
%    dim_u  - dimension of the input
%    n_p    - number of past time steps which are considered
%
% Outputs:
%    nlnARX - generated nonlinearARX object
%
% Example:
%    f = @(y,u) [y(1,1) + u(1,1) - y(2,1); ...
%                   y(3,1) + u(2,1)*cos(y(1,1)); ...
%                   y(5,1) + u(4,1)*sin(y(1,1))];
%    dt = 0.25;
%    sys = nonlinearARX(f,dt,3,2,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT, linearARX

% Authors:       Laura Luetzow
% Written:       24-April-2023
% Last update:   ---
% Last revision: 14-October-2024 (MW, use .n_p instead of .nrOfStates)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    prev_ID = 0;                % previous assigned identifier for polyZon 
    mFile = [];                 % function handle dynamic equation
    jacobian = [];              % function handle jacobian matrix
    hessian = [];               % function handle hessian tensor
    thirdOrderTensor = [];      % function handle third-order tensor
    dt = [];                    % sampling time
    n_p = [];                   % number of past time steps
end

methods
    
    % class constructor
    function nlnARX = nonlinearARX(varargin)
        
        % 0. check number of input arguments
        assertNarginConstructor(2:6,nargin);

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'nonlinearARX')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,fun,dt,dim_y,dim_u,n_p] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,dt,dim_y,dim_u,n_p);   

        % 4. compute properties
        name = char(name);

        % 5a. instantiate parent class
        nlnARX@contDynamics(name,0,dim_u,dim_y);

        % 5b. assign object properties
        nlnARX.dt = dt;
        nlnARX.n_p = n_p;
        nlnARX.mFile = fun;
        nlnARX.jacobian = eval(['@jacobian_',name]);
        nlnARX.hessian = eval(['@hessianTensor_' nlnARX.name]);
        nlnARX.thirdOrderTensor = eval(['@thirdOrderTensor_' nlnARX.name]);

        nlnARX.prev_ID = 10;
    end    

    function setPrevID(nlnARX,id_new)
        nlnARX.prev_ID = id_new;
    end

    function S = generate_indPolyZonotope(nlnARX,S)
        % transform set to a polynomial zonotope with new identifiers

        if ~isa(S,'polyZonotope')
            S = polyZonotope(S);
        end
        id_old = S.id;
        id_new = nlnARX.prev_ID+1 : nlnARX.prev_ID+length(S.id);
        S = replaceId(S, id_old, id_new);

        if ~isempty(id_new)
            setPrevID(nlnARX, id_new(end));
        end
    end


    function R0 = getR0(nlnARX,Y,varargin)
        % stack the previous p output sets to get the NARX 
        % state set for time k+1
        
        narginchk(2,4);
        [type,k] = setDefaultValues({"standard",length(Y)},varargin);

        if isa(Y, 'double')
            yVec = Y;
            Y = cell(size(Y,2),1);
            for j=k-nlnARX.nrOfStates+1:k
                Y{j} = generate_indPolyZonotope(nlnARX, yVec(:,j));
            end
        end

        for j=1:nlnARX.n_p
            if type == "poly" && ~isa(Y{j}, 'polyZonotope')
                Y{j} = generate_indPolyZonotope(nlnARX, Y{j});
            end
            if j==1
                R0 = Y{j};
            else
                if isa(Y{j}, 'polyZonotope')
                    R0 = stack(R0, Y{j});
                else
                    R0 = cartProd(R0, Y{j});
                end
            end
        end

    end

    function [u_stacked, U_stacked] = getStackedU(nlnARX,u,U,varargin)
        % stack the current and previous p inputs sets to get the NARX 
        % input sets for each time point

        narginchk(2,4);
        type = setDefaultValues({"standard"},varargin);

        u_stacked = zeros((nlnARX.n_p+1)*nlnARX.nrOfInputs, size(u,2));
        if nargin == 2
            U_stacked = [];
            for j = nlnARX.n_p+1:size(u,2)
                u_stacked(:,j) = reshape(u(:,j-nlnARX.n_p:j),[],1);
            end
            return
        end

        % compute stacked U-sets 
        if ~iscell(U)
            % U is constant set
            % --> create an indepenent set U for each time step

            U_const = U;
            U = cell(size(u,2),1);
            for j = 1:size(u,2)
                if type == "poly"
                    U{j} = generate_indPolyZonotope(nlnARX, U_const);
                else
                    U{j} = U_const;
                end
            end
        end
        
        U_stacked = cell(size(u,2), 1);
        for j = nlnARX.n_p+1:size(u,2)
            U_stacked_j = U{j} + u(:,j);
            for i=1:nlnARX.n_p
                if isa(U{j-i}, "polyZonotope")
                    U_stacked_j = stack(U{j-i} + u(:,j-i),U_stacked_j);
                else
                    U_stacked_j = cartProd(U{j-i} + u(:,j-i),U_stacked_j);
                end
            end
            U_stacked{j} = U_stacked_j;
            u_stacked(:,j) = center(U_stacked_j);
        end
    end

    % update system dynamics for the new augmented input [u; w] where w is
    % the process noise acting on all states 
    function nlnARX = augment_u_with_w(nlnARX)
        throw(CORAerror('CORA:notSupported',...
            'Not implemented for nonlinearARX.'))
    end

    % update system dynamics for the new augmented input [u; v] where v is
    % the measurement noise acting on all outputs 
    function nlnARX = augment_u_with_v(nlnARX)
        dim_y = nlnARX.nrOfOutputs;
        idz_u_old = [];
        dim_u_old = nlnARX.nrOfInputs;
        dim_u_new = nlnARX.nrOfInputs + dim_y;
        for i = 0:nlnARX.n_p
            idz_u_old = [idz_u_old, i*dim_u_new+1:i*dim_u_new+dim_u_old];
        end
        idz_unew = dim_u_old + 1: dim_u_new;

        nlnARX.mFile = @(y,u) nlnARX.mFile(y,u(idz_u_old)) + u(idz_unew);
        nlnARX.nrOfInputs = dim_u_new;
    end
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [name,fun,dt,dim_y,dim_u,n_p] = aux_parseInputArgs(varargin)

    def_name = 'nonlinearARX';
    fun = @()[]; dt = 0; dim_y = []; dim_u = []; n_p = [];

    % no input arguments
    if nargin == 0
        name = def_name;
        return
    end

    % parse depending on whether first input argument is the name
    if ischar(varargin{1}) || isa(varargin{1},'string')
        [name,fun,dt,dim_y,dim_u,n_p] = ...
            setDefaultValues({def_name,fun,dt,dim_y,dim_u,n_p},varargin);
    else
        name = def_name;
        [fun,dt,dim_y,dim_u,n_p] = ...
            setDefaultValues({fun,dt,dim_y,dim_u,n_p},varargin);
    end
    
end

function aux_checkInputArgs(name,dt,dim_y,dim_u,n_p)

    if CHECKS_ENABLED

        % check name (not empty because default name is not empty)
        % sampling time has to be a scalar larger than zero
        % dim_y and dim_u have to be numeric, scalar integer > 0
        inputArgsCheck({{name,'att',{'char','string'}};...
            {dt,'att','numeric',{'positive','scalar'}};...
            {dim_y,'att','numeric',{'integer','scalar'}};...
            {dim_u,'att','numeric',{'integer','scalar'}};...
            {n_p,'att','numeric','integer','scalar'}});

    end

end

% ------------------------------ END OF CODE ------------------------------
