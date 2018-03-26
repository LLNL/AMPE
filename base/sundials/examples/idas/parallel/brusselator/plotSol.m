% ====================================================================================
% Plot u from Y = (u,v)
% ====================================================================================

function plotSol(x,y,Y)

nx = length(x);
ny = length(y);

[u,v] = Y2UV(nx, ny, Y);

surfc(x,y,u');
shading interp
view(-15,35)
axis tight
box on
grid off
%%xlabel('x','FontSize',16,'FontWeight','bold');
%%ylabel('y','FontSize',16,'FontWeight','bold');

return



function [u,v] = Y2UV(nx,ny,Y)

y2 = reshape(Y, 2, nx*ny);

u = reshape(y2(1,:), ny, nx);
v = reshape(y2(2,:), ny, nx);

return
