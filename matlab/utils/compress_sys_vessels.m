
function [csys, compressinfo] = compress_sys_vessels(sys, nvessmodes, tok)

A = sys.amat;
B = sys.bmat;
nx = size(A,1);

% perform a balanced realization on the vessel currents

% compute balanced realization 
ivess = tok.nc + (1:tok.nv);
Avess = A(ivess,ivess);
Bvess = B(ivess,:);
Cvess = eye(tok.nv, tok.nv);
Pvess = ss(Avess,Bvess,Cvess,0);
[~,~,Tv,~] = balreal(Pvess);
Tvi = inv(Tv);

% format the transformation matrices
iuse = 1:nvessmodes;
Tv = Tv(iuse,:);
Tvi = Tvi(:,iuse);
Tx = blkdiag(eye(tok.nc), Tv);
Txi = blkdiag(eye(tok.nc), Tvi); 

if size(Tx,1) == nx-1    % plasma current is included in states
  Tx = blkdiag(Tx,1);
  Txi = blkdiat(Txi,1);
end


% reduce models
csys.amat = Tx*A*Txi;
csys.bmat = Tx*B;
csys.dpcurrtdx = sys.dpcurrtdx * Txi;
csys.drcurdx = sys.drcurdx * Txi;
csys.dzcurdx = sys.dzcurdx * Txi;

compressinfo = variables2struct(Tv,Tvi,Tx,Txi);





















