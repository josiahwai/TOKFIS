function ts = y2ts(y, iy)

fds = fieldnames(iy);

ts = struct;
for i = 1:length(fds)
  if isnumeric(iy.(fds{i}))    
    ts.(fds{i}) = timeseries(y.Data(:,iy.(fds{i})), y.Time, 'Name', fds{i});
  end
end









