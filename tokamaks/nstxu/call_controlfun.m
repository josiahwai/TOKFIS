
function u = call_controlfun(r, rdot, y, ydot, tnow, CONFIG)

% convert unicode to string. (simulink has shitty string-handling
% capabilities)
controlfun_nam = char(CONFIG.controlfun_unicode);   

% create handle to controlfun
controlfun = str2func(controlfun_nam);

% call controlfun
u = controlfun(r, rdot, y, ydot, tnow, CONFIG);






