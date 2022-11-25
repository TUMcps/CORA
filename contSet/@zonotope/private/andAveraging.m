function [Z] = andAveraging(zonolist,varargin)
% andAveraging - computes the intersection between list of zonotopes
% according to [1]
%
% Syntax:
%    Z = and(zonolist,options)
%
% Inputs:
%    zonolist    - list of zonotopes
%    varargin{1} - method to calculate weights
%                  ('normGen' (default), 'volume', 'radius')
%    varargin{2} - closedform
%                  (true / false)
%    varargin{3} - sumofw (free parameter for the closed form of normGen)
% 
% Outputs:
%    Z - zonotope object enclosing the intersection
%
% Example:
%    zonol{1} = zonotope([2 2 2;1 2 0]);
%    zonol{2} = zonotope([3 1 -1 1;3 1 2 0]);
%    zonol{3} = zonotope([1 3 -1 1;2 3 -2 0]);
%
%    res = andAveraging(zonol);
%
%    figure
%    hold on
%    plot(zonol{1},[1,2],'r');
%    plot(zonol{2},[1,2],'r-+');
%    plot(zonol{3},[1,2],'r-*');
%    plot(res,[1,2],'k');
%    % other supported options for normGen:
%    % find the best sum of w's
%    res = andAveraging(zonol,'normGen',false);
%    % sum of w's =0.9 for the closed form
%    res = andAveraging(zonol,'normGen',true,0.9);
%
%
% References:
%    [1] Amr Alanwar, Jagat Jyoti Rath, Hazem Said, Matthias Althoff
%       Distributed Set-Based Observers Using Diffusion Strategy
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Amr Alanwar
% Written:       9-Feb-2020
% Last update:   9-Mar-2020 (add closed form to normGen)
% Last update:   22-Mar-2020 (add free parameter sum of w's and ability to
%                   switch between closed form and optimization technique)
% Last revision: ---
%
%------------- BEGIN CODE --------------

%2 inputs
if nargin==1 %default
    %The optimization function is based on norm of the generators
    method='normGen';
    sumofw =1;
    closedform = true;
    %3 inputs
elseif nargin==2 % choose other methods
    method =varargin{1};
    if strcmp(method,'normGen')
        sumofw =1;
        closedform = true;
    end
elseif nargin==3 % choose closed form or not for normGen option
    method =varargin{1};
    closedform = varargin{2};
    sumofw =1;
    if ~strcmp(method,'normGen')
        disp('combination is not supported')
        return;
    end
elseif nargin==4 % choose closed form and specific sum of w's
    method =varargin{1};
    closedform = varargin{2};
    if strcmp(method,'normGen') && closedform
        sumofw =varargin{3};
    else
        disp('combination is not supported')
        return;
    end
else
    disp('combination is not supported')
    return;
end


if strcmp(method,'normGen') && closedform
    tVec = zeros(1,length(zonolist));
    w = zeros(1,length(zonolist));
    %find Analytical solution
    invtVecSum = 0;
    for ii=1:length(zonolist)
        tVec(ii)=trace(zonolist{ii}.generators*zonolist{ii}.generators');
        invtVecSum = invtVecSum + 1/tVec(ii);
    end
    
    for ii=1:length(zonolist)
        w(ii)= sumofw/(tVec(ii) * invtVecSum);
    end
elseif strcmp(method,'normGen')|| strcmp(method,'volume') || strcmp(method,'radius')
    w0=(1/length(zonolist))*ones(length(zonolist),1);
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
    %find the weights
    w = fminunc(@fun,w0, options);
else
    disp('option is not supported.');
    return;
end

cen = w(1)*zonolist{1}.center;
gen = w(1)*zonolist{1}.generators;
for ii =2:length(zonolist)
    gen = [gen w(ii)*zonolist{ii}.generators];
    cen = cen + w(ii)*zonolist{ii}.center;
end
cen = cen/sum(w);
gen = (1/sum(w))*gen;
Z = zonotope([cen,gen]);


    function nfro = fun(w)
        gen_inter = w(1)*zonolist{1}.generators;
        cen_inter = w(1)*zonolist{1}.center;
        for i =2:length(zonolist)
            gen_inter = [gen_inter w(i)*zonolist{i}.generators];
            cen_inter = cen_inter + w(i)*zonolist{i}.center;
        end
        gen_inter = (1/sum(w))*gen_inter;
        cen_inter = cen_inter/sum(w);
        
        if strcmp(method,'normGen')
             nfro = norm(gen_inter,'fro');           
        elseif strcmp(method,'radius')
            nfro = radius(zonotope([cen_inter gen_inter]));
        elseif strcmp(method,'volume')
           nfro = volume(zonotope([cen_inter gen_inter]));
        end
        
    end

end
