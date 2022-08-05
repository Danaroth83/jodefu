function g=meanshiftinit(x,h) ;
% initialize the mean-shift procedure

% What is the maximum distance we want to take into account
tolk=1e-3 ; tolx=sqrt(-2*log(tolk)) ;
hmax=tolx * h ;   
g=rangegrid(x,hmax) ;
g.hmax=hmax ;
g.h=h ;

g.plot=0 ;
%g.x=x ;

if g.plot,
  clf ; plot(x(:,2),x(:,1),'b.') ; axis equal ; hold on ;
end ;