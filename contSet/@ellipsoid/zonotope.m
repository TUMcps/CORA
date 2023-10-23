function Z = zonotope(E,varargin)
% zonotope - over-approximates an ellipsoid by a zonotope
%
% Syntax:
%    E = zonotope(E)
%    E = zonotope(E,m)
%    E = zonotope(E,m,mode)
%
% Inputs:
%    E - ellipsoid object
%    m - (optional) number of generators
%    mode - (optional) Specifies whether function uses a lower bound on the 
%           minimum zonotope norm or the exact value:
%           * 'outer:box':      overapprox. parallelotope using
%                               enc_parallelotope
%           * 'outer:norm':     uses enc_zonotope(E,m,'exact') with exact
%                               norm value
%           * 'outer:norm_bnd': not implemented yet (throws error)
%           * 'inner:box':      inner approx. parallelotope using
%                               insc_parallelotope
%           * 'inner:norm'      uses insc_zonotope(E,m,'exact') with exact
%                               norm value
%           * 'inner:norm_bnd': uses insc_zonotope(E,m) with an bound on
%                               the norm value
%           * default:          same as 'outer:norm_bnd'
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Zenc = zonotope(E,10,'outer:norm');
%    Zinsc = zonotope(E,10,'inner:norm');
%    Zbox = zonotope(E);
%
%    figure; hold on;
%    plot(E);
%    plot(Zinsc,[1,2],'r');
%    plot(Zenc,[1,2],'k');
%    plot(Zbox,[1,2],'m');
%
% References:
%    [1] V. Gaßmann, M. Althoff. "Scalable Zonotope-Ellipsoid Conversions
%           using the Euclidean Zonotope Norm", 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       11-October-2019
% Last update:   08-June-2021 (moved handling of degenerate case here)
%                04-July-2022 (VG, class array cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
if isempty(varargin)
    mode = 'outer:box';
    % does not matter, set to any positive integer
    m = 1;
elseif length(varargin)==1
    m = varargin{1};
    % default
    mode = 'outer:norm';
elseif length(varargin)==2
    m = varargin{1};
    mode = varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {mode,'str',{'outer:box','outer:norm','outer:norm_bnd',...
                        'inner:box','inner:norm','inner:norm_bnd'}}});

if isempty(m)
    % only allow if mode is 'outer:box' or 'inner:box'
    if ~strcmp(mode,'inner:box') && ~strcmp(mode,'outer:box')
        throw(CORAerror('CORA:wrongValue','second','integer > 0'));
    end
    % set m to 1 to pass test below
    m = 1;
end

inputArgsCheck({{m,'att',{'numeric'},{'scalar','positive','integer'}}});


% compute rank and dimension of ellipsoid
rankE = rank(E);
dimE = dim(E);
c = center(E);

% handle degenerate case
if ~(rankE == dimE)
    if rankE==0
        Z = zonotope(E.q);
        return;
    end
    % compute svd to find U matrix transforming E.Q to diagonal matrix (to
    % isolate degenerate dimensions)
    [U,Qt,~] = svd(E.Q);
    
    % construct non-degenerate ellipsoid
    E = ellipsoid(Qt(1:rankE,1:rankE));
    
    % construct revert transformation matrix
    T = U(:,1:rankE);
end


switch mode
    case 'outer:box'
        Z = enc_parallelotope(E);
    case 'outer:norm'
        Z = enc_zonotope(E,m,'exact');
    %case 'outer:norm:bnd'% we do not implement the test used to compute the
    %necessary lower bound on Z as in [1] since this test generally does
    %not result in a very good lower bound
    %    Z = enc_zonotope(E,m);
    case 'inner:box'
        Z = insc_parallelotope(E);
    case 'inner:norm'
        Z = insc_zonotope(E,m,'exact');
    case 'inner:norm_bnd'
        Z = insc_zonotope(E,m);
    otherwise
        throw(CORAerror('CORA:wrongValue','third',...
            "be 'outer:box', 'outer:norm', 'inner:box', 'inner:norm', or 'inner:norm_bnd'"));
end

if ~(rankE == dimE)
    Z = T*Z + c;
end

% ------------------------------ END OF CODE ------------------------------
