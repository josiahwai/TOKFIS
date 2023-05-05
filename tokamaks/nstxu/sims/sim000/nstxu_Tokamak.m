function Tokamak = nstxu_Tokamak(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE: Tokamak = nstxu_Tokamak(s)
%
%  PURPOSE: Setup for block Tokamak in nstxu.slx
%
%  INPUTS:  s, simulation settings from setup_nstxu.m
%           tok, TokSys description NSTXU (from connect_nstxu)
%           init, initial equilibrium, used as input to gsinit.m
%
%  OUTPUT: Tokamak, data for subsystem Tokamak in nstxu.slx
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%
%  WRITTEN BY:  Anders Welander ON 2021-01-25
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struct_to_ws(s);


% Begin configuration of the Simulink module GSevolve
% by specifying user-defined options in the variable opt

% HOW TO EVOLVE
% The coil and vessel currents are always evolved according to Vapp = Rs*Is + Mss*Isdot
% Two options exist for evolving the internal profiles
% When evolve_option = 1 there are only three equations that are satisfied for internal profiles,
% therefore the profiles must be constrained to 3 degrees of freedom with constraints = 1
% The three equations are (with input being the left-hand side):
% (1)   Vcd - Rp*Ip = Lp*Ipdot   (current drive voltage - scalar resistance*scalar Ip = inductive plasma voltage)
% (2)   dli/dt      = dli/dt     (d(li)/dt is controlled from the outside with this input)
% (3)   dbetap/dt   = dbetap/dt  (d(betap)/dt is controlled from the outside with this input)
% When evolve_option = 2 the first 2 equations above are replaced with equations for Ohm's parallel law


opt.evolve_option = 1;  % NOTE: shapesim_lite.slx only allows evolve_option==1 
opt.constraints = 1;

% nrepeat can be set to 1 if constraints > 0. It creates corrections to the evolution of Ip (equation 1)
% These can be useful when the evolving plasma is double null. However, this option (and underlying method) may soon be retired.
if opt.constraints > 0
  opt.nrepeat = 0;
end

% OUTPUTS
% The output from the Simulink module GSevolve (called y) is an array with many different signals
% These are specified here in opt.outputs
% A description of the outputs is available in Tokamak.c.iy.info after the setup has finished
opt.outputs = char(...
'gamma',...
'it','ic','iv','sp','sf','rbdef','zbdef',...
'rb','zb','R8','Z8','rsurf','zsurf','aminor','bminor',...
'elong','tril','triu','squo','squi','sqli','sqlo','drsep',...
'psic','psiv','psit',...
'psizr',...
'L1t','jpar','qpsi','pres',...
'ivg',...
'isdiverted',...
'fl','lv','bp','rog',...
'Ip','Ipeq','dpsiplade',...
'vdp','psidp','brdp','bzdp',...
'rcur','zcur',...
'psibry','psimag','rmaxis','zmaxis',...
'Lp','psipla','Atot',...
'cpasma', 'li','betap','Wth','gsphase','fluxerror',...
'now','texec','nupdate','ncalc');

% OTHER SETTINGS
opt.dpsibarlimit = 0.02; % Limit on how much normalized flux can change before new plasma response is calculated
opt.nn = 1;              % number of old response calculations to save in memory (only 1 for now)
opt.dtmin = step_time*0; % Be very careful with using this one, set to zero most of the time
opt.amin = 0;
opt.amax = 0;
opt.rdp = 0.6;   % location at which to measure breakdown voltage
opt.zdp = 0.0;

% Use faeq initiation method
opt.fc.init = 1;
opt.fc.Ip = 100e3;
opt.term.enable = 0;
opt.term.gamma = 10000;
opt.nsavebd = 2;

% the initial R*Bt 
RBtvac = init.rzero * init.bzero;
opt.rzero = init.rzero;         % rzero is an arbitrary reference radius
opt.bzero = init.bzero;         % bzero*rzero is the value that matters, it is RBtvac

if start_time <= 0
  init.fpol = RBtvac;
  init.pres = 0;
end

% Extra resistance (in addition to the PS resistance in MDS and tok.resc), and extra inductance
opt.Rext = 0;
opt.Lext = 0;

% Create the configuration data for GSevolve
Tokamak.c = gsconfig(tok,opt);
Tokamak.c.plots = gsevolve_plot_opts(tok);

% e = gsinit(Tokamak.c, init);
% [Tokamak.d, Tokamak.e, r, b] = gsupdate(e.x, e, Tokamak.c, 1);
% Tokamak.yinit = Tokamak.d.y0 + Tokamak.d.C*e.x;
% Tokamak.yinit(isnan(Tokamak.yinit)) = 0;

[Tokamak.d, Tokamak.e, r, b] = gsupdate(init.x, init, Tokamak.c, 1);
Tokamak.yinit = Tokamak.d.y0 + Tokamak.d.C*init.x;
Tokamak.yinit(isnan(Tokamak.yinit)) = 0;


%--------------------------------------------------------------
% particles
Tokamak.particles.vbd = 1;
Tokamak.particles.tbd = 0e-3;
Tokamak.particles.iy = Tokamak.c.iy;


%--------------------------------------------------------------
% pressure

opt.info.evolve_option = Tokamak.c.info.evolve_option;
opt.evolve_option = Tokamak.c.evolve_option;

opt.info.ndpdt = 'Dimension of signal dpdt';
opt.ndpdt = numel(Tokamak.c.iu.sp);

opt.info.P = 'Proportional gain for correcting the error betap-ref.betap';
opt.P = 100;

opt.info.sp0 = 'Tokamak.c.sp0';
opt.sp0 = Tokamak.c.sp0;

ref.info.betap = 'Preprogrammed poloidal beta';
ref.betap = betap;

Tokamak.pressure.iy = Tokamak.c.iy;
Tokamak.pressure.opt = opt;
Tokamak.pressure.ref = ref;


%--------------------------------------------------------------
% current

% Initial plasma
Tokamak.c.faeq.li = max(0.4, li.Data(1));
Tokamak.c.faeq.betap = max(0.01, betap.Data(1));

opt.info.evolve_option = Tokamak.c.info.evolve_option;
opt.evolve_option = Tokamak.c.evolve_option;

opt.info.nVp = 'Dimension of signal Vp';
opt.nVp = numel(Tokamak.c.iu.sf);

opt.info.P = 'Proportional gain for correcting the error li-ref.li';
opt.P = 100;


% Preprogrammed values will be in ref
ref = [];

ref.info.Rp = 'preprogrammed scalar resistance';
ref.Rp = Rp;

ref.info.li = 'Preprogrammed li';
ref.li = li;

Tokamak.current.iy = Tokamak.c.iy;
Tokamak.current.opt = opt;
Tokamak.current.ref = ref;


%----------------------------------------------------------------
% documentation

Tokamak.info.c = 'Output from gsconfig, used to configure GSevolve with tok info plus options';
Tokamak.info.initsim = 'Data from simulation between start of voltage and start of plasma';
Tokamak.info.d = 'Output from gsupdate with dynamics objects such as ABCD for start_time';
Tokamak.info.e = 'Output from gsupdate with equilibrium info for start_time';
Tokamak.info.yinit = 'Output from GSevolve at start_time';
Tokamak.info.particles = 'Setup data for subsystem particles';
Tokamak.info.pressure = 'Setup data for subsystem pressure';
Tokamak.info.current = 'Setup data for subsystem current';

