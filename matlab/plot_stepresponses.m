% Helper function for plotting data stored as a struct of timeseries (or
% timeseries-like, with Time and Data fields). 


function fighandle = plot_nyquists(G, idx, titles, nrows, fighandle, varargin)

N = length(idx);
if ~exist('nrows','var') || isempty(nrows), nrows = min(2,N); end
if ~exist('fighandle','var') || isempty(fighandle), fighandle = figure; end
if isempty(varargin), varargin = {}; end

ncols = ceil(N/nrows);

figure(fighandle);

for i = 1:N
  subplot(nrows, ncols, i)
  hold on
  k = idx(i);
  nyquist(G(k,k), varargin{:});
  title(titles{i}, 'fontsize', 16)
  axis([-4 1 -2 2])
  drawnow
end

  






























