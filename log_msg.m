function log_msg(level, msg)
%LOG_MSG Simple logging utility
%   log_msg('error','something bad');
%   log_msg('info','normal info');
%
%   Level is optional; default 'info'. Messages are printed with timestamp.
    if nargin < 2
        msg = level;
        level = 'info';
    end
    if ~ischar(level)
        level = 'info';
    end
    t = datestr(now,'yyyy-mm-dd HH:MM:SS');
    fprintf('[%s] %s: %s\n', t, upper(level), msg);
end
