function out = num2sepstr(numin, format, sep)
% NUM2SEPSTR Convert to string with separation at thousands.
%
% out = NUM2SEPSTR(numin,[format],[sep]) formats numin to a string
%   according to the specified format ('%f' by default) and adds the
%   sepcified thousands seperators (commas by default).
%
% For non-scalar numin, num2sepstr outpts a cell array of the same shape
%   as numin where num2sepstr is called on each value in numin.
%
% String length from format, when specified, is applied before commas are
%   added; Instead of...
%
%   >> num2sepstr(1e6,'% 20.2f') % length = 22
%   ans =
%       '          1,000,000.00'
% 
% ...try...
% 
%   >> sprintf('% 20s',num2sepstr(1e6,'%.2f')) % length = 20
%   ans =
%       '        1,000,000.00'
%
% See also SPRINTF, NUM2STR
%
% Created by:
%   Robert Perrotta

% Thanks to MathWorks community members Stephen Cobeldick and Andreas J.
% for suggesting a faster, cleaner implementation using regexp and
% regexprep.

if nargin < 2
    format = ''; % we choose a format below when we know numin is scalar and real
end
if nargin < 3
    sep = ',';
end

if numel(numin)>1
    out = cell(size(numin));
    for ii = 1:numel(numin)
        out{ii} = num2sepstr(numin(ii), format, sep);
    end
    return
end

if ~isreal(numin)
    out = sprintf('%s+%si', ...
        num2sepstr(real(numin), format, sep), ...
        num2sepstr(imag(numin), format, sep));
    return
end

autoformat = isempty(format);
if autoformat
    if isinteger(numin) || mod(round(numin, 4), 1) == 0
        format = '%.0f';
    else
        format = '%.4f'; % 4 digits is the num2str default
    end
end

str = sprintf(format, numin);

if isempty(str)
    error('num2sepstr:invalidFormat', ...
        'Invalid format (sprintf could not use "%s").', format)
end

out = regexpi(str, '^(\D*\d{0,3})(\d{3})*(\D\d*)?$', 'tokens', 'once');
if numel(out)
    out = [out{1}, regexprep(out{2}, '(\d{3})', [sep,'$1']), out{3}];
else
    out = str;
end

if autoformat
    % Trim trailing zeros after the decimal. (By checking above for numbers
    % that look like integers using autoformat, we avoid ever having ONLY
    % zeros after the decimal. There will always be at least one nonzero
    % digit following the decimal.)
    out = regexprep(out, '(\.\d*[1-9])(0*)', '$1');
end

end
