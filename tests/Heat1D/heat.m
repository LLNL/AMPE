function heat

ntintervals=5
tspan=5000.

global L;
L = 1.;
global nx;
nx = 60;
global alpha;
alpha=1.e-5

%specifies slab symmetry
m = 0;

%time intervals where solution desired
x = linspace(0.,L,nx);
t = linspace(0,tspan,ntintervals+1)

sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);

temperature = sol(:,:,1); %temperature

fht = figure(1);
surf(x,t,temperature);
title('T(x,t)');
xlabel('Distance x');
ylabel('Time t');
saveas(fht, 'temperature.png', 'png');

timeslots = linspace(1,ntintervals,ntintervals)

colors=['k','r','g','b','c','m']

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

global alpha;

s=[0.];
f=[alpha*DuDx(1)];
c=[1.];

% --------------------------------------------------------------------------
%initial conditions
function u0 = pdex4ic(x)

Tinit=2.*x+sin(2.*pi*x)+1.;

u0 = [Tinit];                                 

% --------------------------------------------------------------------------
%BC: p(x,t,u)+q(x,t)*flux=0
%phase u(1) is 1 on left and 0 on right
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
pl = [ul(1)-1.];                               
ql = [0. ];                                  
pr = [-2.];                            
qr = [ 1.e5];                                  
