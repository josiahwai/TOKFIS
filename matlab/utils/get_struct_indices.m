function iy = get_struct_indices(s, fds)

iy = struct;

idx = 0;
for i = 1:length(fds)
  fd = fds{i};
  n = size(s.(fd), 1);  
  idx = idx(end)+1:idx(end)+n;  
  iy.(fd) = idx;
end 
