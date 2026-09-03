[x,y] = meshgrid(linspace(-2,2,100));
z = x.^2 + y.^2;

surf(x,y,z)
[F,V] = surf2patch(x,y,z,'triangles');
TR = triangulation(F,V);

stlwrite(TR,'paraboloid.stl');