function plot_coils(tok, c, alpha)

nc = size(tok.fcdata,2);
fcdata = tok.fcdata;

% housekeeping
if ~exist('c','var') || isempty(c)
  c = [1 1 1] * 0.75;
end

if ischar(c)
  c = bitget(find('krgybmcw'==c)-1,1:3);
end

if ~exist('alpha','var') || isempty(alpha)
  alpha = 1;
end

if size(c,1) ~= nc, c = ones(nc,1)*c; end
if isscalar(alpha), alpha = alpha*ones(nc,1); end


% plot coils
for i = 1:nc
  
  y = fcdata(1,i) - fcdata(3,i)/2;
  x = fcdata(2,i) - fcdata(4,i)/2;
  dx = fcdata(4,i);
  dy = fcdata(3,i);
  position = [x y dx dy];

  h = rectangle('Position', position);
  h.FaceColor = [c(i,:) alpha(i)];
  h.EdgeColor = c(i,:) * min(1-alpha(i), 0.8);

end


% axis
xmax = max(fcdata(2,:));
xmin = min(fcdata(2,:));
ymax = max(fcdata(1,:));
ymin = min(fcdata(1,:));
axis equal
dx = (xmax-xmin)*0.05;
dy = (ymax-ymin)*0.05;

axis([xmin-dx xmax+dx ymin-dy ymax+dy])





