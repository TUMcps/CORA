function writeSparseMatrixOptimized(ind,var,fid,tayMod)

    % loop over all non-empty entries
    for i = 1:length(ind.row)
        if ~tayMod
            str = [var,'(',num2str(ind.row(i)),',',num2str(ind.col(i)), ...
                   ') = out(',num2str(ind.index(i)),');'];
        else
            str = [var,'(',num2str(ind.row(i)),',',num2str(ind.col(i)), ...
                   ') = interval(out(',num2str(ind.index(i)),'));'];
        end
          
       fprintf(fid, '%s\n', str);
    end
end