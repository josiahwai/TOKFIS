% Use LQR to design the coil current controller. This can often be 
% more effective at decoupling the control than PID

function K = nstxu_coilcontrol_lqr(tok, sys)


% to have "balanced" coil responses, weights are proportional to 
% circuit inductance squared. Could also tune weights individually
L = diag(tok.mcc);   
wt.ic = (L*180).^2;        % coil weights
wt.iv = zeros(tok.nv,1);   % vessel weights are zero


% LQR weighting matrices
R = eye(tok.nc);  
Q = diag([wt.ic; wt.iv]);  


% solve LQR
K = lqr(sys.amat, sys.bmat, Q, R);


% K has feedback terms on both coil and vessel currents. However, vessel
% currents are usually not measured. We will just delete the vessel current
% part and see how controller performs
K = K(:,1:tok.nc);



