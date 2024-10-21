function xrank = ndfsort_rank1(x)
%NDFSORT_RANK1 Fast non dominated sorting algorithm for 2 objectives only - returns only rank 1 solutions.
%
%   XRANK = NDFSORT_RANK1(X) performed on a (N x 2) matrix X where N is the
%   number of solutions and there are 2 objectives separates out the
%   solutions into different levels (ranks) of non-domination. The returned
%   vector XRANK contains only rank 1 solutions. Solutions with rank 1 come
%   from the non-dominated front of X.
%
%   Remarks: 
%
%   Based on the sorting method described on page 184 of: "A fast and
%   elitist multiobjective genetic algorithm: NSGA-II" by Kalyanmoy Deb,
%   Amrit Pratap, Sameer Agarwal, T. Meyarivan. IEEE Transactions on
%   Evolutionary Computation Vol. 6, No. 2, April 2002, pp. 182-197.
%
%   (c) Dominic Searson 2009-2015
%
%   GPTIPS 2
%
%   See also SELECTION, GPMODELFILTER, PARETOREPORT, POPBROWSER

P = size(x,1);
np = zeros(P,1); %current domination level
xrank = np;  %rank vector

for p=1:P
    
    for q=1:P
        
        if q ~= p
            if doesx1Domx2(x(q,:),x(p,:))
                np(p) = np(p) + 1;
            end
        end
        
    end
    
    if np(p) == 0
        xrank(p) = 1;
    end
    
end

function result = doesx1Domx2(x1,x2)
%Returns true if first solution dominates second, false otherwise. Here,
%dominance means if one solution beats the other on 1 criterion and is no
%worse on the other criterion (not strict dominance).

result = false;

if (x1(1) < x2(1) && x1(2) <= x2(2)) || (x1(1) <= x2(1) && x1(2) < x2(2))
    result = true;
end
