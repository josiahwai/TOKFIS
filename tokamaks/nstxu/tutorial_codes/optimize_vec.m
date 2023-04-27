% Given the plasma current distribution and info about target boundary,
% estimate the coil currents/applied flux.
%
% Goal is to minimize:
%   J = (ytarg - y)' * W * (ytarg - y) 
%  where y is the measured variables (field at locations, flux, etc), 
%  ytarg is the target measurements, W is a weighting matrix. 
% 
% Given that y = G*I where I is the coil and vessel currents and G is the
% greens functions from current sources to measurements. 
% 
% And subject to constraints on I, such as current limits etc. This problem
% is rewritten as min J = I'*G'*W*G*I - 2*ytarg'*W*I subject to the
% constraints and can be solved with quadprog. 
%

% ===============
% Define geometry
% ===============

function vec = optimize_vec(tok, targ, wt, opts)

% Housekeeping...

% optimization defaults
default.ub = inf(tok.nc,1);
default.lb = -inf(tok.nc,1);
default.Aineq = [];
default.bineq = [];
default.Aeq = [];
default.beq = [];
if ~exist('opts','var'), opts = struct; end
opts = copyfields(default, opts, [], 1);
struct_to_ws(opts);

% defaults for targ and wt
targdefault.ic  = zeros(tok.nc,1);
targdefault.br  = zeros(size(targ.r));
targdefault.bz  = zeros(size(targ.r));
targdefault.psi = zeros(size(targ.r));

wtdefault.ic  = ones(size(targdefault.ic));
wtdefault.br  = ones(size(targdefault.br));
wtdefault.bz  = ones(size(targdefault.bz));
wtdefault.psi = ones(size(targdefault.psi));

targ = copyfields(targdefault, targ, [], 1);
wt   = copyfields(wtdefault, wt, [], 1);


% Now on to the real optimization...

% compute responses
C.ic = eye(tok.nc);
C.br = gridresponse2pt(tok.rg, tok.zg, tok.gbr2c, targ.r, targ.z);
C.bz = gridresponse2pt(tok.rg, tok.zg, tok.gbz2c, targ.r, targ.z);
C.psi = gridresponse2pt(tok.rg, tok.zg, tok.mpc, targ.r, targ.z);

fds = {'ic', 'br', 'bz', 'psi'};

G = struct2vec(C, fds);
ytarg = struct2vec(targ, fds);
W = diag(struct2vec(wt, fds));

f = -G' * W * ytarg;
H = G'*W*G;
H = (H+H')/2;


% solve quadprog
vec = quadprog(H, f, Aineq, bineq, Aeq, beq, lb, ub);





