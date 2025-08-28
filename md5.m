function h = md5(str)
%MD5 Compute MD5 hash of input string, returning lowercase hex digest.
%   h = MD5(str) returns the lowercase MD5 hash for the given string.
    if isstring(str)
        str = char(str);
    end
    h = lower(hash('md5', str, 'hex'));
end
