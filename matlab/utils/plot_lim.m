

function plot_lim(tok)


if size(tok.limdata,1) ~= 2
  tok.limdata = tok.limdata';
end

hold on
plot(tok.limdata(2,:), tok.limdata(1,:), 'color', [1 1 1] * 0.4, ...
  'linewidth', 1.5);

axis equal
axis([min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)])








