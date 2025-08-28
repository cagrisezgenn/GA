function badVars = validate_finite(varNames)
%VALIDATE_FINITE Check that variables are finite and non-empty.
%   BADVARS = VALIDATE_FINITE(VARNAMES) takes a cell array of variable
%   names (as strings) evaluated in the caller workspace. It returns a cell
%   array BADVARS containing the names of any variables that are empty or
%   contain non-finite values (NaN or Inf).

badVars = {};
for k = 1:numel(varNames)
    name = varNames{k};
    try
        val = evalin('caller', name);
    catch
        badVars{end+1} = name; %#ok<AGROW>
        continue;
    end
    if isempty(val) || any(~isfinite(val(:)))
        badVars{end+1} = name; %#ok<AGROW>
    end
end
end
