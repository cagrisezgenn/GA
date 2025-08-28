function val = psa_cache(action, key, val_in)
%PSA_CACHE Persistent cache for PSA results.
%   val = PSA_CACHE(action, key, val_in)
%   Actions:
%     'get'  - retrieve value for key (empty if not found)
%     'set'  - store value and return it
%     'clear'- clear all cached values
    persistent cache
    if isempty(cache)
        cache = containers.Map('KeyType','char','ValueType','any');
    end
    switch action
        case 'get'
            if isKey(cache, key)
                val = cache(key);
            else
                val = [];
            end
        case 'set'
            cache(key) = val_in;
            val = val_in;
        case 'clear'
            cache = containers.Map('KeyType','char','ValueType','any');
            val = [];
        otherwise
            error('psa_cache:unknownAction', 'Unknown action %s', action);
    end
end
