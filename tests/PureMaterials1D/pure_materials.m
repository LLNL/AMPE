%1D test problem for Beckermann et al. model 
%(Beckermann et al., J. Comput. Phys. 124, 468-496 (1999) )
%units: m,K,s,J
function pure_materials

ntintervals=4
tspan=40.

global L;
L = 1.8e-3;
global nx;
nx = 600;
global delta;
delta=1.27e-5/6.;

%specifies slab symmetry
m = 0;

x = linspace(0.,L,nx);

%time intervals where solution desired
t = linspace(0,tspan,ntintervals+1)

global Tm;
Tm = 933.6
global ml;
ml=-260.;

sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);

phi = sol(:,:,1); %phase
temperature = sol(:,:,2); %temperature

fh = figure(1);
surf(x,t,phi);
title('phi(x,t)');
xlabel('Distance x');
ylabel('Time t');
saveas(fh, 'phi.png', 'png');

fht = figure(2);
surf(x,t,temperature);
title('T(x,t)');
xlabel('Distance x');
ylabel('Time t');
saveas(fht, 'temperature.png', 'png');

timeslots = linspace(1,ntintervals,ntintervals)

colors=['k','r','g','b','c','m']

fhp = figure(4);

for timeslot = 1:ntintervals
  p=phi(timeslots(timeslot),:);
  plot(x,p,colors(timeslot))
  hold on

  %save data in text files
  a=[x',p'];
  fid = fopen(sprintf('phase%d.dat', timeslot),'w');
  if fid<0
    disp('ERROR in opening file!!!')
    exit
  end
  fprintf(fid,'%10.7f  %8.3f\n',a');
  fclose(fid);
end
hold off

saveas(fhp, 'profile.png', 'png');

fht = figure(5);

for timeslot = 1:ntintervals
  p=temperature(timeslots(timeslot),:);
  plot(x,p,colors(timeslot))
  hold on

  %save data in text files
  a=[x',p'];
  fid = fopen(sprintf('temperature%d.dat', timeslot),'w');
  if fid<0
    disp('ERROR in opening file!!!')
    exit
  end
  fprintf(fid,'%4.2f  %8.3f\n',a');
  fclose(fid);
end
hold off

saveas(fht, 'profileT.png', 'png');


% --------------------------------------------------------------------------
%define fluxes and sources for all three coupled equations
function [c,f,s] = pdex4pde(x,t,u,DuDx)

phi=u(1);
phis=min(max(0,phi),1);

global Tm;
global ml;
global delta;
muk=2.6e-5;
Gamma=2.41e-7;
K=9.5e8/2.58e6; %L/cP
alpha=100.e-6; %Al: 97 mm2/s, Copper: 111 mm2/s

%mT=(Tm-u(2)+ml*Cl);
mT=(Tm-u(2)+ml*0.048);

d2=delta*delta;
phisrc=muk*( -Gamma*phi*(1.-phi)*(1.-2.*phi)/d2+mT*phi*(1.-phi)/delta);

s=[phisrc;
   K*phisrc];

fs=muk*Gamma*DuDx(1);
f=[fs;
   alpha*DuDx(2)+K*fs];

c=[1.;1.];

% --------------------------------------------------------------------------
%initial conditions
function u0 = pdex4ic(x)
global Tm;
global nx;
global L;
global ml;

h=L/nx;
invdelta=1./h;
d=x-0.2*L;
v = 0.5*(1.+tanh(-0.5*d*invdelta));

cl=0.048;
Tinit=Tm+cl*ml;

u0 = [v; Tinit];                                 

% --------------------------------------------------------------------------
%BC: p(x,t,u)+q(x,t)*flux=0
%phase u(1) is 1 on left and 0 on right
%zero flux for u(3) (composition)
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
pl = [ul(1)-1.; -1.e-2];                               
ql = [0.;        1.;  ];                                  
pr = [ur(1);     1.e-4];                            
qr = [0.;        1.;  ];                                  
