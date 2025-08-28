function key = psa_hash(t, ag, T_vec, zeta)
%PSA_HASH Compute a hash string for PSA input combinations.
%   key = PSA_HASH(t, ag, T_vec, zeta) returns an MD5 hash based on the
%   provided vectors/scalars. This is used to cache PSA computations.
    data = [];
    if ~isempty(t),     data = [data; t(:)]; end
    if ~isempty(ag),    data = [data; ag(:)]; end
    if ~isempty(T_vec), data = [data; T_vec(:)]; end
    if ~isempty(zeta),  data = [data; zeta(:)]; end
    key = hash('md5', char(typecast(data,'uint8')).');
end
