function res = testLongDuration_quadZonotope_reduceRedistGirard()
% testLongDuration_quadZonotope_reduceRedistGirard - unit test for zonotope 
%                                                    reduction with the 
%                                                    method "redistGirard"
%
% Syntax:  
%    res = testLongDuration_quadZonotope_reduceRedistGirard()
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

% RANDOM TESTS

for d = 3:5
   for e = 5:10
       for g = 6:10
          for j = 1:3
              % select option
              if j == 1
                 options.reduceRedistGirard = 'ascend'; 
              elseif j == 2
                 options.reduceRedistGirard = 'descend';
              else
                 options.reduceRedistGirard = 'none';
              end                  

              % construct random zonotope
              E = rand(d,e)-ones(d,e);
              G = rand(d,g)-ones(d,g);
              c = rand(d,1);
              zono = zonotope([c,E,G]);
              vert = vertices(zono);

              % reduce zonotope
              n = e*(e-1)/2;
              qZ = quadZonotope(c,E,zeros(d,e),zeros(d,n),G);
              qZ = reduce(qZ,'redistGirard',1,options);
              c_ = center(qZ);
              [G_,~,~,Grest_] = generators(qZ);
              zonoRed = zonotope([c_,G_,Grest_]);

              % check if the reduction is overapproximative
              [C,D] = halfspaceRepresentation(zonoRed);

              for i = 1:size(vert,2)
                  temp = ((C*vert(:,i)-D) <= 1e-10);
                  if ~all(temp);
                      res = false;
                      disp('test_quadZonotope_reduceGirard failed');
                      return;
                  end
              end
          end
       end
    end
end
    
end


function [C,d] = halfspaceRepresentation(zono)

    %obtain number of generators, dimensions
    %Z=deleteAligned(Z);
    Z=get(zono,'Z');
    c=Z(:,1);
    G=Z(:,2:end);
    [dim,nrOfGenerators]=size(G);

    if dim > 1
        %get number of possible facets
        comb=combinator(nrOfGenerators,dim-1,'c');

        %build C matrices
        C=[];
        for i=1:length(comb(:,1))
            indices=comb(i,:);
            Q=G(:,indices);
            v=ndimCross(Q);
            C(end+1,:)=v'/norm(v);
        end

        %remove NaN rows due to rank deficiency
        index = find(sum(isnan(C),2));
        if ~isempty(index)
            C(index,:) = [];
        end
    else
        C = G;
    end

    %build d vector
    %determine delta d
    deltaD=zeros(length(C(:,1)),1);
    for iGen=1:nrOfGenerators
        deltaD=deltaD+abs(C*G(:,iGen));
    end

    %compute dPos, dNeg
    dPos=C*c+deltaD;
    dNeg=-C*c+deltaD;

    C = [C;-C];
    d = [dPos;dNeg];
    
end