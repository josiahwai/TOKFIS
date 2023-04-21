


function circ = sparc_circ()

% coil to circuit mappings
cccirc = [1 1 2 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 -19];
Pcc = cccirc_to_Pcc(cccirc);
nc = max(abs(cccirc));

% circuit names
ccnames = {'CS1U', 'CS1L', 'CS2U', 'CS2L', 'CS3U', 'CS3L', 'PF1U', ...
  'PF1L', 'PF2U', 'PF2L', 'PF3U', 'PF3L', 'PF4U', 'PF4L', 'DV1U', ...
  'DV1L', 'DV2U', 'DV2L', 'VS1'}';


% power supply voltage limits
vmax = [2.6 2.6 1 1 1 1 1 1 1 1 2.6 2.6 2.6 2.6 0.55 0.55 0.55 0.55 1.1]' * 1e3; 
vmin = -vmax;

icmax = [45*ones(1,14) 32*ones(1,4) 10]' * 1e3;
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















































































