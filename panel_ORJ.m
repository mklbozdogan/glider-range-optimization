function [xm,vtang,Cp,Cl,Cd,Cm] = panel_ORJ(x,y,alpha,xmom)

np = numel(x) ;
nodtot = np-1 ;

alpi = alpha * pi()/ 180.0 ;

% --- Compute slopes and arc length
xm = ( x(2:end) + x(1:end-1) ) / 2.0 ;
ym = ( y(2:end) + y(1:end-1) ) / 2.0 ;
dx = diff(x) ;
dy = diff(y) ;
dist = sqrt(dx.^2 + dy.^2) ;
sinthe = dy ./ dist ;
costhe = dx ./ dist ;

% --- Compute a steady solution at AOA alpi
pi2inv = 0.5 / pi() ;

cosalf = cos(alpi) ;
sinalf = sin(alpi) ;

kutta = nodtot + 1 ;
A = zeros(kutta,kutta) ;
b = zeros(kutta,1) ;
aan = zeros(nodtot,nodtot) ;
bbn = zeros(nodtot,nodtot) ;

for i = 1:nodtot
for j = 1:nodtot
  if j == i
    flog = 0.0 ;
    ftan = pi ;
  else
    dxj = xm(i) - x(j) ;
    dxjp = xm(i) - x(j+1) ;
    dyj = ym(i) - y(j) ;
    dyjp = ym(i) - y(j+1) ;
    flog = log((dxjp^2 + dyjp^2) / (dxj^2 + dyj^2)) / 2.0 ;
    ftan = atan2(dyjp*dxj - dxjp*dyj, dxjp*dxj + dyjp*dyj) ;
  end

  ctimtj = costhe(i)*costhe(j) + sinthe(i)*sinthe(j) ;
  stimtj = sinthe(i)*costhe(j) - costhe(i)*sinthe(j) ;
              
  A(i,j) = pi2inv * (ftan*ctimtj + flog*stimtj) ;
  bbn(i,j) = pi2inv * (flog*ctimtj - ftan*stimtj) ;
  A(i, kutta) = A(i, kutta) + bbn(i,j) ;

  if (i == 1) || (i == nodtot)
    A(kutta, j) = A(kutta, j) - bbn(i,j) ;
    A(kutta, kutta) = A(kutta, kutta) + A(i,j) ;
  end
  aan(i,j) = A(i,j) ;
end
b(i) = sinthe(i)*cosalf - costhe(i)*sinalf ;
end
b(kutta) = -(costhe(1) + costhe(nodtot))*cosalf ...
           -(sinthe(1) + sinthe(nodtot))*sinalf ;

q = A\b ;
gamma = q(kutta) ;
q = q(1:nodtot) ;
       
vtang = cosalf*costhe + sinalf*sinthe - bbn*q + gamma*aan*ones(nodtot,1) ;
Cp = 1.0 - vtang.^2 ;

cfx =  Cp'*dy ;
cfy = -Cp'*dx ; 

Cl = cfy*cosalf - cfx*sinalf ;
Cd = cfx*cosalf + cfy*sinalf ;
Cm = Cp' * ( dx.*(xm-xmom) + dy.*ym ) ;

