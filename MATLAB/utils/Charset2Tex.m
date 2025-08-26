function str = Charset2Tex(Charset)
str = cell(size(Charset));

for i = 1:length(Charset)
    
    str{i} = strrep(Charset{i}, 'x', 'x_');
    
end