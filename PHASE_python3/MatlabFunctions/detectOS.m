function os = detectOS()

% detectOS Determines the operating system.
%
% Output:
%   os: A string representing the operating system:
%       - 'windows' 
%       - 'linux'
%       - 'macos_intel'
%       - 'macos_apple_silicon' 

% Check if it's Windows
if ispc
    os = 'windows';
    return;
end

% Check if it's macOS
if ismac
    [~, output] = system('sysctl -n machdep.cpu.brand_string'); 
    if contains(output, 'Apple')
        os = 'macos_apple_silicon';
    else
        os = 'macos_intel';
    end
    return;
end

% Check if it's Linux
if isunix && ~ismac
    os = 'linux'; 
end
end