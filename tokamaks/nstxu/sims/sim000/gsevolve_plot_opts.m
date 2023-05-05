% Some default options for plots used in gsevolve sims
% Anything that is specified in opts will overwrite the defaults

function plots = gsevolve_plot_opts(tok, opts)


% PLOT SETTINGS
plots.updates = 0;
plots.dtplot = 0.002;
plots.tplot = [-1:0.001:0.2 0.21:0.005:9999];
plots.axis = [min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)];

% Boundary
plots.boundary.LineWidth = 2;
plots.boundary.Color = [0 0.9 0];
plots.boundary.patch = false;
plots.boundary.EdgeColor = [.8 .6 .1];
plots.boundary.FaceColor = [1 1 .8];

% Boundary-defining point
plots.bdef.LineWidth = 0;
plots.bdef.Color = [1 .4 0];
plots.bdef.Marker = 'x';
plots.bdef.MarkerSize = 18;

% Touch point
plots.touch.LineWidth = 3;
plots.touch.Color = [1 .8 0];
plots.touch.Marker = '*';
plots.touch.MarkerSize = 12;

% Diagnostic points
plots.dp.LineWidth = 0;
plots.dp.Color = [0 1 0];
plots.dp.Marker = '.';
plots.dp.MarkerSize = 10;

% Flux contours
plots.flux_contours.LineWidth = 1;
plots.flux_contours.nvac = 0;
plots.flux_contours.npla = 0;
plots.flux_contours.psibar = linspace(0,2,9);

% Camera
plots.camera.showit = 0; % 1 = show camera image of plasma
plots.camera.animate = 1; % 0 = show still of initial plasma
plots.camera.nxpix = 64;
plots.camera.nzpix = 128;

% Limiter
if plots.camera.showit
  plots.limiter.LineWidth = 6;
else
  plots.limiter.LineWidth = 3;
end
if plots.camera.showit
  plots.limiter.Color = [0 1 0];
else
  plots.limiter.Color = [0.5, 0.5, 0.5];
end
plots.limiter.patch = false;
plots.limiter.EdgeColor = [0.5 0.5 0.5];
plots.limiter.FaceColor = [0.7 0.7 0.7];

% Vessel
plots.vessel.LineWidth = 2;
plots.vessel.Color = [.1 .1 .8];
plots.vessel.patch = false;
plots.vessel.EdgeColor = [0 1 1];
plots.vessel.FaceColor = [0 0 1];

% Coils
plots.coils.LineWidth = 1;
plots.coils.Color = [1 .1 0];
plots.coils.patch = false;
plots.coils.EdgeColor = [0 1 1];
plots.coils.FaceColor = [0 0 1];

% Text
plots.text.FontSize = 18;
plots.save_frames = 0;

if exist('opts','var') && isstruct(opts)
  plots = copyfields(plots, opts, [], 1);
end




