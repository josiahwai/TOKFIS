
function circ = nstxu_circ()

% coil to circuit mappings
cccirc = [1 2 0 0 3 3 4 4 4 4 0 0 0 0 0 0 5 5 5 5 6 6 6 6 7 7 0 0 8];
Pcc = cccirc_to_Pcc(cccirc);
nc = max(abs(cccirc));


% circuit names
ccnames = {'OH', 'PF1AU', 'PF2U', 'PF3U', 'PF5', 'PF3L', 'PF2L', 'PF1AL'}';


% power supply voltage limits
vmax = [4048 1012 2024 2024 3036 2024 2024 1012]';
vmin = -vmax;


% power supply current limits
icmax = [20  15 15 8 0 8 15 15]' * 1e3;
icmin = [-20 0 0 -13 -24 -13 0 0]' * 1e3;

% the iy struct holds indices of each circuit so that they can be
% referred to easily by name 
for i = 1:nc
  iy.(ccnames{i}) = i;
end


% descriptions
d.Pcc = 'mapping from coils to circuits: icx=Pcc*ic, ic=pinv(Pcc)*icx';
d.ccnames = 'circuit names';
d.nc = 'number of circuits';
d.cccirc = 'circuit connection vector';
d.vmax = 'power supply max voltage';
d.vmin = 'power supply min voltage';
d.icmax = 'max circuit current';
d.icmin = 'min circuit current';
d.iy = 'circuit names and indices';
descriptions = d;


% save circuit data
circ = variables2struct(cccirc, Pcc, nc, ccnames, vmax, vmin, icmax, ...
  icmin, iy);


% sort circuit data
fds               = sort(fields(circ));
circ              = reorderstructure(circ, fds{:});
descriptions      = reorderstructure(descriptions, fds{:});
circ.descriptions = descriptions; 















































































