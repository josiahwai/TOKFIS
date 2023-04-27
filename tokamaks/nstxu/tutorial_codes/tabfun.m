% similar in spirit to cellfun, apply a function to every tab in a tabgroup

function tabfun(figh, fun, varargin)


for i = 1:length(figh.Children)
  h = figh.Children(i);
  if isa(h, 'matlab.ui.container.TabGroup')
    tabgroup = h;
    break
  end
end


tabs = tabgroup.Children;

for i = 1:length(tabs)
  ax = tabs(i).Children(1);
  axes(ax);
  fun(varargin{:})
end



















