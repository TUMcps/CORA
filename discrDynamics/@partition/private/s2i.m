function indices=s2i(Obj,subscripts)
    multiplicator = Obj.nrOfSegments';
    for i = 2:length(multiplicator)
        multiplicator(i) = multiplicator(i)*multiplicator(i-1);
    end
    multiplicator = [1, multiplicator(1:end-1)];
    
    indices = (multiplicator * (subscripts - 1)')+1;
    indices=indices.*(prod((subscripts<=Obj.nrOfSegments').*(subscripts>0),2))';
    

%     for iDim=1:length(currentIndex)
%         sizeMX=size(MX,1) ; %<-- don't think you really need to do this but oh well
%         MX=[repmat(MX,length(currentIndex{iDim}),1),reshape(ones(sizeMX,1)*(currentIndex{iDim}),[],1)];
%     end
%     Multiples=ones(length(currentIndex),1);
%     for i = 1:(length(currentIndex)-1)
%         Multiples(i,1)=prod(Obj.nrOfSegments((i+1):end));
%     end
%     segments=(MX-1)*Multiples;
end