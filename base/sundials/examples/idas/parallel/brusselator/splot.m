nx = 41;
ny = 41;


Nx = 2*nx;
Ny = 2*ny;


x = linspace(0.0,1.0,Nx);
y = linspace(0.0,1.0,Ny);

%% Local grids

x0 = x(1:nx);
y0 = y(1:ny);

x1 = x(nx+1:end);
y1 = y(1:ny);

x2 = x(1:nx);
y2 = y(ny+1:end);

x3 = x(nx+1:end);
y3 = y(ny+1:end);

%% Initial conditions

Y0 = load('ic_0a.txt');
Y1 = load('ic_1a.txt');
Y2 = load('ic_2a.txt');
Y3 = load('ic_3a.txt');

figure;

set(gcf,'position',[600 300 600 600])
hold on

plotSol(x0,y0,Y0);
plotSol(x1,y1,Y1);
plotSol(x2,y2,Y2);
plotSol(x3,y3,Y3);

%%title('u(0.0) \equiv u_0','FontSize',16);

set(gca,'FontName','times')
%%set(gca,'FontWeight','bold')
set(gca,'FontSize',20)


%% Final solution

Y0 = load('sol_0a.txt');
Y1 = load('sol_1a.txt');
Y2 = load('sol_2a.txt');
Y3 = load('sol_3a.txt');

figure;

set(gcf,'position',[700 400 600 600])
hold on

plotSol(x0,y0,Y0);
plotSol(x1,y1,Y1);
plotSol(x2,y2,Y2);
plotSol(x3,y3,Y3);

%%title('u(1.0)','FontSize',16);

set(gca,'FontName','times')
%%set(gca,'FontWeight','bold')
set(gca,'FontSize',20)



%% adjoint variables

Y0 = load('adj_0a.txt');
Y1 = load('adj_1a.txt');
Y2 = load('adj_2a.txt');
Y3 = load('adj_3a.txt');

figure;

set(gcf,'position',[800 500 600 600])
hold on

plotSol(x0,y0,Y0);
plotSol(x1,y1,Y1);
plotSol(x2,y2,Y2);
plotSol(x3,y3,Y3);

%%title('dg(t_f)/d_u_0 \equiv \lambda(0.0)','FontSize',16);


set(gca,'FontName','times')
%%set(gca,'FontWeight','bold')
set(gca,'FontSize',20)


