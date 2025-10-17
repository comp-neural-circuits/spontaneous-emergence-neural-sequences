function  arr=retrieve_sp(sparsearr)
        arr = zeros(sparsearr.m, sparsearr.n);
        for j = 1:sparsearr.n
                index = (sparsearr.colptr(j) : (sparsearr.colptr(j+1) - 1));
                arr(sparsearr.rowval(index), j) = sparsearr.nzval(index);
        end
end