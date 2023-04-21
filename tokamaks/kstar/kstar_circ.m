
function circ = kstar_circ()

% coil to circuit mappings
cccirc = [1 2 3 4 5 6 7 1 2 8 9 10 11 7 13 12 -13 12];
Pcc = cccirc_to_Pcc(cccirc);
nc = max(abs(cccirc));


% circuit names
ccnames = {'PF1UL', 'PF2UL', 'PF3U', 'PF4U', 'PF5U', 'PF6U', 'PF7UL', ...
  'PF3L', 'PF4L', 'PF5L', 'PF6L', 'IRCUL', 'IVCUL'}';


% voltage limits
vmax = [1 1 0.5 0.5 1 1 1 0.5 0.5 1 1 0.5 1]' * 1e3;
vmin = -vmax;


% current limits
icmax = [7.5 8.5 9 9 7 3.5 3.5 9 9 7 3.5 8 3]' * 1e3; 
icmin = -icmax;


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



































