% Helper function for plotting data stored as a struct of timeseries (or
% timeseries-like, with Time and Data fields). 


function fighandle = plot_step_responses(G, idx, titles, tmax, nrows, fighandle, varargin)

N = length(idx);
if ~exist('nrows','var') || isempty(nrows), nrows = min(2,N); end
if ~exist('fighandle','var') || isempty(fighandle), fighandle = figure; end
if isempty(varargin), varargin = {}; end

ncols = ceil(N/nrows);

figure(fighandle);

for i = 1:N
  subplot(nrows, ncols, i)
  k = idx(i);
  step(G(k,k), tmax)
  hold on
  title(titles{i}, 'fontsize', 16)
  drawnow
end

  






























