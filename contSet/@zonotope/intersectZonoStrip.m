function Zres = intersectZonoStrip(z1,hl,Rl,yl,varargin)
% intersectZonoStrip - computes the intersection between one zonotope and
%    list of strips according to [1]
%    the strip is defined as | hx-y | <= d
%
% Syntax:  
%    Zres = intersectZonoStrip(z1,hl,Rl,yl,varargin)
%
% Inputs:
%    z1 - zonotope object
%    h1 - 
%    R1 - 
%    y1 -
%    varargin - methods to calculate the weights
%               'normGen' default and has analytical solution
%               'svd'
%               'radius'
%
% Outputs:
%    res - boolean whether obj is contained in Z, or not
%
% Example: (three strips and one zonotope)
%    hl{1} = [1 0];
%    Rl{1} = 5;
%    yl{1} = -2;
% 
%    hl{2} = [0 1];
%    Rl{2} = 3;
%    yl{2} = 2;
% 
%    hl{3} = [1 1];
%    Rl{3} = 3;
%    yl{3} = 2;
% 
%    Z = zonotope([1 2 2 2 6 2 8;1 2 2 0 5 0 6 ]);
%    res_zono = intersectZonoStrip(Z,hl,Rl,yl);
% 
%    % just for comparison
%    poly = mptPolytope([1 0;-1 0; 0 1;0 -1; 1 1;-1 -1],[3;7;5;1;5;1]);
%    Zpoly = Z & poly;
% 
%    figure; hold on 
%    plot(Z,[1 2],'r-+');
%    plot(poly,[1 2],'r-*');
%    plot(Zpoly,[1 2],'b-+');
%    plot(res_zono,[1 2],'b-*');
% 
%    legend('zonotope','strips','zono&poly','zonoStrips');
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
% Written:       09-Mar-2020
% Last update:   ---              
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin==4
    %The optimization function is based on norm of the generators
    method='normGen';
elseif nargin==5
    method =varargin{1};
end

    H = generators(z1);
if strcmp(method,'svd') || strcmp(method,'radius') 
    lambda0=zeros(length(z1.center),length(Rl));
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
    %find the weights
    lambda = fminunc(@fun,lambda0, options);
elseif strcmp(method,'normGen')
    % Find the analytical solution
    h_combined=[];
    for i=1:length(hl)
        h_combined = [ h_combined ; hl{i}];
    end    
    gamma=eye(length(hl));
    num= H*H'*h_combined';
    den = h_combined * H*H' * h_combined' ;
    for i=1:length(hl)
        den = den + gamma(:,i) *Rl{i}^2* gamma(:,i)';
    end
    
    lambda = num * den^-1;
else
    disp('Method is not supported');
    return;
end



%prepare center
c_new=z1.center;
for i=1:length( Rl)
    c_new = c_new + lambda(:,i)*( yl{i} - hl{i}*z1.center );
end

%prepare generators
part1 = eye(length(z1.center));
for ii=1:length(Rl)
    part1 = part1 - lambda(:,ii)*hl{ii};
    part2(:,ii) = Rl{ii}*lambda(:,ii);
end
part1 = part1 * H;
H_new = [part1 part2];
Zres = zonotope([c_new H_new]);



    function nfro = fun(lambda)
        part1 = eye(length(z1.center));
        for ii=1:length(Rl)
            part1 = part1 - lambda(:,ii)*hl{ii};
            part2(:,ii) = Rl{ii}*lambda(:,ii);
        end
        part1 = part1 * H;
        H_new = [part1 part2];
        if strcmp(method,'svd')
            nfro = sum(svd(H_new));
        elseif strcmp(method,'radius')
            nfro = radius(zonotope([zeros(length(z1.center),1) H_new]));
        end
        
    end


end
