function v = getfield_default(S, fname, defaultVal)
%GETFIELD_DEFAULT Safe field access with default value.
%   v = GETFIELD_DEFAULT(S, fname, defaultVal) returns the value of the field
%   FNAME from struct S if it exists and is non-empty. If the field is
%   missing, empty, or contains non-finite numeric values, DEFAULTVAL is
%   returned instead.
%
%   This helper mirrors the more robust implementation originally embedded
%   in viscous.m and is now available as a standalone utility for reuse.

    if ~isstruct(S) || ~isfield(S, fname) || isempty(S.(fname))
        v = defaultVal;
        return;
    end
    val = S.(fname);
    if isnumeric(val)
        if isempty(val) || any(~isfinite(val(:)))
            v = defaultVal;
        else
            v = val;
        end
    else
        % For non-numeric types, only emptiness is checked.
        v = val;
    end
end
