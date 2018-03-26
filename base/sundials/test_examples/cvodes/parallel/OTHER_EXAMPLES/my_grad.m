clear;

grad000;
grad001;
grad002;
grad003;


figure;
hold on
trans = 0.7;
ecol  = 'none';
xp=[4.000000 16.000000];
yp=[8.000000 12.000000];
%zp=[8.000000 12.000000];
zp=[];
ns = length(xp)*length(yp)*length(zp);


[X,Y,Z]=meshgrid(x0,y0,z0);
s0=slice(X,Y,Z,g0,xp,yp,zp);
for i = 1:length(s0)
  set(s0(i),'FaceAlpha',trans);
  set(s0(i),'EdgeColor',ecol);
end


[X,Y,Z]=meshgrid(x1,y1,z1);
s1=slice(X,Y,Z,g1,xp,yp,zp);
for i = 1:length(s1)
  set(s1(i),'FaceAlpha',trans);
  set(s1(i),'EdgeColor',ecol);
end


[X,Y,Z]=meshgrid(x2,y2,z2);
s2=slice(X,Y,Z,g2,xp,yp,zp);
for i = 1:length(s2)
  set(s2(i),'FaceAlpha',trans);
  set(s2(i),'EdgeColor',ecol);
end


[X,Y,Z]=meshgrid(x3,y3,z3);
s3=slice(X,Y,Z,g3,xp,yp,zp);
for i = 1:length(s3)
  set(s3(i),'FaceAlpha',trans);
  set(s3(i),'EdgeColor',ecol);
end
view(3)

shading interp
axis equal



figure;
hold on
iso = 1.0;
trans = 0.3;
fcol = 'yellow';

[X,Y,Z]=meshgrid(x0,y0,z0);
p = patch(isosurface(X,Y,Z,g0,iso));
isonormals(X,Y,Z,g0,p);
set(p, 'FaceColor', fcol, 'EdgeColor', 'none', 'FaceAlpha', trans);
daspect([1 1 1])

[X,Y,Z]=meshgrid(x1,y1,z1);
p = patch(isosurface(X,Y,Z,g1,iso));
isonormals(X,Y,Z,g1,p);
set(p, 'FaceColor', fcol, 'EdgeColor', 'none', 'FaceAlpha', trans);
daspect([1 1 1])


[X,Y,Z]=meshgrid(x2,y2,z2);
p = patch(isosurface(X,Y,Z,g2,iso));
isonormals(X,Y,Z,g2,p);
set(p, 'FaceColor', fcol, 'EdgeColor', 'none', 'FaceAlpha', trans);
daspect([1 1 1])


[X,Y,Z]=meshgrid(x3,y3,z3);
p = patch(isosurface(X,Y,Z,g3,iso));
isonormals(X,Y,Z,g3,p);
set(p, 'FaceColor', fcol, 'EdgeColor', 'none', 'FaceAlpha', trans);
daspect([1 1 1])


view(3)
camlight; lighting phong


