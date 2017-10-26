%1D test problem for Beckermann et al. model 
%(Beckermann et al., J. Comput. Phys. 124, 468-496 (1999) )
%units: m,K,s,J
function pure_materials

ntintervals=5
tspan=50.
topJperum3=1.e-6;

global L;
L = 400.;
global nx;
nx = 256;
nx=512;
global Vm;
Vm=1.e-5
global delta;
delta=1.27e1/6.;
muk=2.6e1;
Gamma=2.41e-1;
global LA;
LA=9.5e3
LA=LA/Vm;
LA=LA*topJperum3
global Tmc;
Tmc=921.12

global deltaT;
deltaT =0.25

theoretical_velocity = muk*deltaT

global M;
global epsilon;
global omega;
M=muk*Tmc/(6.*delta*LA)
epsilon=sqrt(6.*delta*LA*Gamma/Tmc)
omega=((epsilon/delta)^2 )/32.

%specifies slab symmetry
m = 0;

x = linspace(0.,L,nx);

%time intervals where solution desired
t = linspace(0,tspan,ntintervals+1)
dt=tspan/ntintervals

%options=odeset('AbsTol',1e-4)
%sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,options);
sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);

phi = sol(:,:,1); %phase

fh = figure(1);
surf(x,t,phi);
title('phi(x,t)');
xlabel('Distance x');
ylabel('Time t');
saveas(fh, 'phi.png', 'png');

timeslots = linspace(1,ntintervals+1,ntintervals+1)

colors=['k','r','g','b','c','m']
fraction=0.05*L;
fhp = figure(4);

for timeslot = 1:ntintervals+1
  p=phi(timeslots(timeslot),:);
  new_fraction=trapz(x,p);
  velocity=(new_fraction-fraction)/dt
  fraction=new_fraction;

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


% --------------------------------------------------------------------------
%define fluxes and sources for all three coupled equations
function [c,f,s] = pdex4pde(x,t,u,DuDx)

global Tmc;
global M;
global epsilon;
global omega;
global LA;
global deltaT;

phi=u(1);
phisrc=-M*32.*omega*phi*(1.-phi)*(1.-2.*phi)+M*6.*LA*deltaT*phi*(1.-phi)/Tmc;

s=[phisrc];

fs=M*epsilon*epsilon*DuDx(1);
f=[fs];

c=[1.];

% --------------------------------------------------------------------------
%initial conditions
function u0 = pdex4ic(x)
global L;
global delta;

invdelta=1./delta;
d=x-0.05*L;
v = 0.5*(1.+tanh(-0.5*d*invdelta));

u0 = [v];                                 

% --------------------------------------------------------------------------
%BC: p(x,t,u)+q(x,t)*flux=0
%phase u(1) is 1 on left and 0 on right
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
pl = [ul(1)-1.];                               
ql = [0.;     ];                                  
pr = [ur(1);  ];                            
qr = [0.;     ];                                  
