
s.vacuum     = 0;
s.start_time = 0;
s.stop_time  = 0.95;
s.step_time  = 5e-6;
s.dt2save    = 50e-6;
s.decimation = round(s.dt2save/s.step_time);
s.init       = init;
s.tok        = tok;

t = linspace(0,1);
s.betap = timeseries(init.betap, t);
s.li    = timeseries(init.li, t);
s.Rp    = timeseries(5e-7, t);

% clear gsupdate cache
clear gsupdate

Tokamak = nstxu_Tokamak(s);
Actuators = nstxu_Actuators();
c = Tokamak.c;
e = Tokamak.e;

