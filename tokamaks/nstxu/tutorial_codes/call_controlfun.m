
function upad = call_controlfun(r, rdot, y, ydot, tnow, CONFIG, padlen)

% read the name of the controlfun, which was passed in as a vector of
% unicode values (not a normal string, b/c simulink has shitty bugs)
controlfun_nam = char(CONFIG.controlfun_unicode);   

% create handle to controlfun
controlfun = str2func(controlfun_nam);

% call controlfun
u = controlfun(r, rdot, y, ydot, tnow, CONFIG);

% add zero-padding so that the output is a fixed length
upad = zeros(padlen,1);
upad(1:CONFIG.nu) = u;




