% match a cellstr of ccnames to indices

function idx = ccnames2idx(ccnames, tok)

tok.ccnames = cellstr(tok.ccnames);
ccnames = cellstr(ccnames);

if any(strcmp(ccnames, 'all'))
  ccnames = tok.ccnames;
end

idx = [];
for i = 1:length(ccnames)
  idx(end+1) = find(strcmpi(ccnames{i}, tok.ccnames));
end

