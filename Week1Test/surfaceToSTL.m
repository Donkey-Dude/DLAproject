function surfaceToSTL(filename, bottomFun, topFun, Rmax, nr, nt)
% surfaceToSTL
%
% Creates a closed 3D-printable STL between two surfaces:
%
%       z = bottomFun(x,y)
%       z = topFun(x,y)
%
% over a circular domain:
%
%       x^2 + y^2 <= Rmax^2
%
% INPUTS:
%   filename  - output STL filename
%   bottomFun - lower surface function @(x,y)
%   topFun    - upper surface function @(x,y)
%   Rmax      - radius of domain
%   nr        - radial resolution
%   nt        - angular resolution
%
% EXAMPLE:
%   surfaceToSTL('object.stl', ...
%                @(x,y) x.^2 + y.^2, ...
%                @(x,y) 16 + 0*x, ...
%                4, 80, 160);

%% -----------------------------
% Generate vertices
% ------------------------------

V = [];

% ----- Bottom center -----
x = 0;
y = 0;

V(end+1,:) = [x, y, bottomFun(x,y)];
bottomCenter = size(V,1);

% ----- Bottom rings -----
bottomRings = zeros(nr,nt);

for i = 1:nr

    r = Rmax * i/nr;

    for j = 1:nt

        % Do NOT include theta = 2*pi separately
        theta = 2*pi*(j-1)/nt;

        x = r*cos(theta);
        y = r*sin(theta);
        z = bottomFun(x,y);

        V(end+1,:) = [x,y,z];
        bottomRings(i,j) = size(V,1);

    end
end


% ----- Top center -----
x = 0;
y = 0;

V(end+1,:) = [x, y, topFun(x,y)];
topCenter = size(V,1);

% ----- Top rings -----
topRings = zeros(nr,nt);

for i = 1:nr

    r = Rmax * i/nr;

    for j = 1:nt

        theta = 2*pi*(j-1)/nt;

        x = r*cos(theta);
        y = r*sin(theta);
        z = topFun(x,y);

        V(end+1,:) = [x,y,z];
        topRings(i,j) = size(V,1);

    end
end


%% -----------------------------
% Generate triangular faces
% -----------------------------

F = [];

%% Bottom center fan

for j = 1:nt

    jnext = mod(j,nt) + 1;

    F(end+1,:) = [
        bottomCenter, ...
        bottomRings(1,jnext), ...
        bottomRings(1,j)
    ];

end


%% Bottom surface

for i = 1:nr-1

    for j = 1:nt

        jnext = mod(j,nt) + 1;

        a = bottomRings(i,j);
        b = bottomRings(i,jnext);
        c = bottomRings(i+1,j);
        d = bottomRings(i+1,jnext);

        F(end+1,:) = [a,d,c];
        F(end+1,:) = [a,b,d];

    end
end


%% Top center fan

for j = 1:nt

    jnext = mod(j,nt) + 1;

    F(end+1,:) = [
        topCenter, ...
        topRings(1,j), ...
        topRings(1,jnext)
    ];

end


%% Top surface

for i = 1:nr-1

    for j = 1:nt

        jnext = mod(j,nt) + 1;

        a = topRings(i,j);
        b = topRings(i,jnext);
        c = topRings(i+1,j);
        d = topRings(i+1,jnext);

        F(end+1,:) = [a,c,d];
        F(end+1,:) = [a,d,b];

    end
end


%% -----------------------------
% Connect outside boundary
% -----------------------------

for j = 1:nt

    jnext = mod(j,nt) + 1;

    b1 = bottomRings(nr,j);
    b2 = bottomRings(nr,jnext);

    t1 = topRings(nr,j);
    t2 = topRings(nr,jnext);

    F(end+1,:) = [b1,b2,t2];
    F(end+1,:) = [b1,t2,t1];

end


%% -----------------------------
% Create triangulation
% -----------------------------

TR = triangulation(F,V);


%% -----------------------------
% Preview
% -----------------------------

figure
trisurf(TR, ...
    'FaceColor',[0.7 0.8 1], ...
    'EdgeColor','none');

axis equal
xlabel('x')
ylabel('y')
zlabel('z')

camlight
lighting gouraud

title('STL Preview')


%% -----------------------------
% Export
% -----------------------------

stlwrite(TR,filename);

fprintf('Created: %s\n',filename);

end
