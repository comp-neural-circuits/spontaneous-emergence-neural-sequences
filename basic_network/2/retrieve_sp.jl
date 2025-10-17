using SparseArrays

function retrieve_sp(sparsearr::Any)
    arr = zeros(sparsearr.m, sparsearr.n)
    for j = 1:sparsearr.n
        index = (sparsearr.colptr[j]:(sparsearr.colptr[j+1]-1))
        arr[sparsearr.rowval[index], j] = sparsearr.nzval[index]
    end
    return arr
end

function retrieve_sp(sparsearr::Dict)
    arr = zeros(sparsearr["m"], sparsearr["n"])
    for j = 1:sparsearr["n"]
        index = (sparsearr["colptr"][j]:(sparsearr["colptr"][j+1]-1))
        arr[sparsearr["rowval"][index], j] = sparsearr["nzval"][index]
    end
    return arr
end
