function result = configsEqual(left,right)
%CONFIGSEQUAL Compare Model configurations including NaN and datetime values.

result = isequaln(orderfields(left),orderfields(right));
end
