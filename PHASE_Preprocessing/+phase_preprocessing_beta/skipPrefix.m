function prefix = skipPrefix()
%SKIPPREFIX Comment out a skipped generated step on each shell platform.

if ispc
    prefix = 'rem ';
else
    prefix = '# ';
end
end
