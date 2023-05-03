function ts = y2ts(y, iy)

fds = fieldnames(iy);

ydata = squeeze(y.Data);
if size(ydata,1) ~= length(y.Time)
  y.Data = ydata';
end


ts = struct;
for i = 1:length(fds)
  if isnumeric(iy.(fds{i}))    
    ts.(fds{i}) = timeseries(y.Data(:,iy.(fds{i})), y.Time, 'Name', fds{i});
  end
end









