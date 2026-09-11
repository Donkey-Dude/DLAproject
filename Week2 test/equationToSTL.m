function equationToSTL(filename, inside, bounds, resolution)

% Create coordinate vectors
x = linspace(bounds(1,1), bounds(1,2), resolution);
y = linspace(bounds(2,1), bounds(2,2), resolution);
z = linspace(bounds(3,1), bounds(3,2), resolution);

% Create 3D grid
[X,Y,Z] = meshgrid(x,y,z);

% Evaluate whether each point is inside the object
solid = inside(X,Y,Z);

% Convert logical values to numbers for isosurface
solid = double(solid);

% Extract boundary between inside (1) and outside (0)
[F,V] = isosurface(X,Y,Z,solid,0.5);

% Check that a surface was actually found
if isempty(F)
    error("No surface found. Check equations, bounds, and resolution.");
end

% Create triangulation
TR = triangulation(F,V);

% Preview
figure
trisurf(TR)

axis equal
xlabel("x")
ylabel("y")
zlabel("z")
title("STL Preview")

% Export
stlwrite(TR,filename);

% Confirm export
if isfile(filename)
    fprintf("STL successfully created:\n%s\n", ...
        fullfile(pwd,filename));
else
    error("STL export failed.");
end

end