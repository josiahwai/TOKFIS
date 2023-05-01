
function upad = call_controlfun(r, rdot, y, ydot, tnow, CONFIG)

% convert unicode to string 
controlfun_nam = char(CONFIG.controlfun_unicode);   

% create handle to controlfun
controlfun = str2func(controlfun_nam);

% 
% u = zeros(CONFIG.nu,1);
% coder.varsize('u', CONFIG.nu, 0);

% call controlfun
u = controlfun(r, rdot, y, ydot, tnow, CONFIG);

upad = zeros(100,1);
upad(1:CONFIG.nu) = u;




