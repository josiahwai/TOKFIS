
function plot_eq(eq, tok, varargin)

try
  rg = eq.rg;
  zg = eq.zg;
catch
  rg = tok.rg;
  zg = tok.zg;
end
plot_lim(tok)
hold on
contour(rg, zg, eq.psizr,eq.psibry * [1 1], varargin{:});












