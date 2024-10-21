function Z = priv_andAveraging(Z_cell,varargin)
% priv_andAveraging - computes the intersection between list of zonotopes
%    according to [1]
%
% Syntax:
%    Z = priv_andAveraging(Z_cell)
%    Z = priv_andAveraging(Z_cell,method)
%    Z = priv_andAveraging(Z_cell,method,closedform)
%    Z = priv_andAveraging(Z_cell,method,closedform,sumofw)
%
% Inputs:
%    Z_cell - cell-array of zonotopes
%    method - method to calculate weights:
%             'normGen' (default)
%             'volume'
%             'radius'
%    closedform - true/false
%    sumofw - free parameter for the closed form of normGen
% 
% Outputs:
%    Z - zonotope object enclosing the intersection
%
% Example:
%    Z1 = zonotope([2 2 2;1 2 0]);
%    Z2 = zonotope([3 1 -1 1;3 1 2 0]);
%    Z = and(Z1,Z2,'averaging');
%
%    figure; hold on;
%    plot(Z1); plot(Z2);
%    plot(Z,[1,2],'k');
%    % % other supported options for normGen:
%    % % find the best sum of w's
%    % res = andAveraging(zonol,'normGen',false);
%    % % sum of w's =0.9 for the closed form
%    % res = andAveraging(zonol,'normGen',true,0.9);
%
% References:
%    [1] Amr Alanwar, Jagat Jyoti Rath, Hazem Said, Matthias Althoff.
%        Distributed Set-Based Observers Using Diffusion Strategy
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Amr Alanwar
% Written:       09-February-2020
% Last update:   09-March-2020 (add closed form to normGen)
% Last update:   22-March-2020 (sum of w's and switch between closed form/optimization)
% Last revision: 06-October-2024 (MW, complete refactor)

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,4);
[method,closedform,sumofw] = setDefaultValues({'normGen',true,1},varargin);
inputArgsCheck({{Z_cell,'att','cell'}; ...
                {method,'str',{'normGen','radius','volume'}}; ...
                {closedform,'att','logical'}; ...
                {sumofw,'att','numeric',{'scalar','nonnegative'}}});

% nargin == 2:
% ...if method == 'normGen', set closedform = true, sumofw = 1, otherwise nothing
% nargin == 3:
% ...method must be 'normGen'
% nargin == 4:
% ...method must be 'normGen' and closedform must be true

% compute weighting factors
if strcmp(method,'normGen') && closedform
    % analytical solution
    tVec = cellfun(@(Z) trace(Z.G*Z.G'),Z_cell,'UniformOutput',true);
    inv_tVecSum = sum(1 ./ tVec);
    w = sumofw ./ (tVec * inv_tVecSum);

else
    % find the weights via optimization
    options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');
    w0 = ones(numel(Z_cell),1)/numel(Z_cell);
    w = fminunc(@aux_fun,w0,options);
end

[c,G] = aux_catWeighted(Z_cell,w);
Z = zonotope(c,G);

    % nested function so that aux_fun has only 1 input argument for fminunc
    function nfro = aux_fun(w)
        
    % add all centers and concatenate all generators from zonotopes
    [cen_inter,gen_inter] = aux_catWeighted(Z_cell,w);
    
    if strcmp(method,'normGen')
        nfro = norm(gen_inter,'fro');           
    elseif strcmp(method,'radius')
        nfro = radius(zonotope(cen_inter, gen_inter));
    elseif strcmp(method,'volume')
        nfro = volume(zonotope(cen_inter, gen_inter));
    end
        
    end

end


% Auxiliary functions -----------------------------------------------------

function [c,G] = aux_catWeighted(zonolist,w)
% add all centers and concatenate all generator matrices from all zonotopes
% in the cell-array, including a weighting factor

n = dim(zonolist{1});
c = zeros(n,1); G = zeros(n,0);
for i = 1:length(zonolist)
    G = [G, w(i)*zonolist{i}.generators];
    c = c + w(i)*zonolist{i}.center;
end
c = c/sum(w);
G = (1/sum(w))*G;

end

% ------------------------------ END OF CODE ------------------------------
