function x = h5read_sparsematrix(file, var)
        x = struct();
        x.m = h5read(file, ['/', var, '/m']);
        x.n = h5read(file, ['/', var, '/n']);
        x.colptr = h5read(file, ['/', var, '/colptr']);
        x.rowval = h5read(file, ['/', var, '/rowval']);
        x.nzval = h5read(file, ['/', var, '/nzval']);
end