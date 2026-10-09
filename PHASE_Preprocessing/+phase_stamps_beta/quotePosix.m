function value = quotePosix(inputValue)
%QUOTEPOSIX Quote one argument for the /bin/bash commands used on Unix.
single = char(39);
embedded = [single char(34) single char(34) single];
value = [single strrep(char(string(inputValue)),single,embedded) single];
end
