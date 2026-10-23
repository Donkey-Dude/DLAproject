function [paths, distance, info] = cubeShortestPaths(p, doPlot, tol)
% [paths,d,info] = cubeShortestPaths([x y z],doPlot,tol)


if nargin >= 1 && (ischar(p) || isa(p,'string'))
    if strcmpi(p,'test'), selftest(); paths={}; distance=[]; info=[]; return; end
    error('Text input must be ''test''.');
end
if nargin < 1 || isempty(p), p=[0.5 1 0]; end
if nargin < 2 || isempty(doPlot), doPlot=true; end
if nargin < 3 || isempty(tol), tol=1e-10; end
validateattributes(tol,{'numeric'},{'real','finite','scalar','positive','<',1e-3});
validateattributes(p,{'numeric'},{'real','finite','vector','numel',3});
p=double(p(:).'); s=[0.5 0 1];
if any(p < -tol | p > 1+tol) || ~any(abs(p)<tol | abs(p-1)<tol)
    error('p must lie on the surface of [0,1]^3.');
end
p(abs(p)<tol)=0; p(abs(p-1)<tol)=1;

V=[0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
F=[1 2 3 4;5 6 7 8;1 2 6 5;4 3 7 8;1 4 8 5;2 3 7 6];
adj=false(6); for a=1:6, for b=1:6
    adj(a,b)=numel(intersect(F(a,:),F(b,:)))==2;
end, end
sf=incident(s); pf=incident(p); sequences={};
for a=sf, visit(a); end
candidates={}; lens=[]; labels={};
for k=1:numel(sequences)
    seq=sequences{k}; n=numel(seq); maps=cell(1,n);
    M=nan(8,2); M(F(seq(1),:),:)=[0 0;1 0;1 1;0 1]; maps{1}=M;
    for j=2:n
        old=F(seq(j-1),:); new=F(seq(j),:); h=intersect(old,new);
        A=maps{j-1}(h(1),:); B=maps{j-1}(h(2),:); e=B-A;
        out=[-e(2) e(1)]; c=mean(maps{j-1}(old,:),1);
        if dot(out,c-A)>0, out=-out; end
        M=nan(8,2); M(h,:)=[A;B];
        for v=setdiff(new,h)
            along=dot(V(v,:)-V(h(1),:),V(h(2),:)-V(h(1),:));
            M(v,:)=A+along*e+out;
        end
        maps{j}=M;
    end
    S=project(s,seq(1),maps{1}); P=project(p,seq(n),maps{n}); D=P-S;
    ts=zeros(1,n+1); ts(end)=1; valid=true;
    for j=1:n-1
        h=intersect(F(seq(j),:),F(seq(j+1),:));
        [lo,hi,ok]=hingeInterval(S,D,maps{j}(h(1),:),maps{j}(h(2),:),tol);
        ts(j+1)=max(ts(j),lo);
        if ~ok || ts(j+1)>hi+tol, valid=false; break; end
        ts(j+1)=min(1,max(0,ts(j+1)));
    end
    if ~valid, continue; end
   
    Q=zeros(n+1,3);
    for j=1:n
        f=F(seq(j),:); m=maps{j}(f,:);
        for t=ts(j:j+1)
            q=S+t*D;
            if any(q<min(m,[],1)-tol | q>max(m,[],1)+tol), valid=false; end
        end
        Q(j,:)=unproject(S+ts(j)*D,seq(j),maps{j});
    end
    if ~valid, continue; end
    Q(abs(Q)<tol)=0; Q(abs(Q-1)<tol)=1;
    Q(end,:)=p; Q(1,:)=s; Q=simplify(Q,tol);
    candidates{end+1}=Q; lens(end+1)=norm(D); labels{end+1}=seq; %#ok<AGROW>
end
if isempty(lens), error('No valid development found.'); end
distance=min(lens); paths={}; faceSequences={};
for k=find(abs(lens-distance)<=tol)
    Q=candidates{k}; duplicate=false;
    for j=1:numel(paths)
        if isequal(size(Q),size(paths{j})) && max(abs(Q(:)-reshape(paths{j},[],1)))<=10*tol
            duplicate=true; break;
        end
    end
    if ~duplicate, paths{end+1}=Q; faceSequences{end+1}=labels{k}; end %#ok<AGROW>
end
info=struct('source',s,'destination',p,'tolerance',tol, ...
    'faceNames',{{'z=0','z=1','y=0','y=1','x=0','x=1'}}, ...
    'faceSequences',{faceSequences},'sequencesExamined',numel(sequences), ...
    'validDevelopments',numel(lens),'pathLengths',cellfun(@curveLength,paths));
if doPlot
    figure('Color','w'); patch('Vertices',V,'Faces',F,'FaceColor',[.8 .85 .9], ...
        'FaceAlpha',.13,'EdgeColor',[.3 .3 .3]); hold on;
    colors=lines(numel(paths));
    for j=1:numel(paths)
        Q=paths{j}; plot3(Q(:,1),Q(:,2),Q(:,3),'-o','LineWidth',2,'Color',colors(j,:));
    end
    plot3(s(1),s(2),s(3),'ko','MarkerFaceColor','g','MarkerSize',9);
    plot3(p(1),p(2),p(3),'ks','MarkerFaceColor','r','MarkerSize',9);
    text(s(1),s(2),s(3),'  s'); text(p(1),p(2),p(3),'  p');
    axis equal; axis([0 1 0 1 0 1]); grid on; view(3);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('%d distinct shortest paths; length %.12g',numel(paths),distance));
    rotate3d on;
end
    function fs=incident(q)
        fs=find([q(3)==0,q(3)==1,q(2)==0,q(2)==1,q(1)==0,q(1)==1]);
    end
    function visit(seq)
        if any(seq(end)==pf), sequences{end+1}=seq; end
     
        for b=find(adj(seq(end),:))
            if ~any(seq==b), visit([seq b]); end
        end
    end
    function q=project(x,f,M)
        ids=F(f,:); o=V(ids(1),:); u=V(ids(2),:)-o; v=V(ids(4),:)-o;
        q=M(ids(1),:)+dot(x-o,u)*(M(ids(2),:)-M(ids(1),:)) ...
            +dot(x-o,v)*(M(ids(4),:)-M(ids(1),:));
    end
    function x=unproject(q,f,M)
        ids=F(f,:); o=M(ids(1),:); u=M(ids(2),:)-o; v=M(ids(4),:)-o;
        x=V(ids(1),:)+dot(q-o,u)*(V(ids(2),:)-V(ids(1),:)) ...
            +dot(q-o,v)*(V(ids(4),:)-V(ids(1),:));
    end
end

function [lo,hi,ok]=hingeInterval(S,D,A,B,tol)
axisAlong=find(abs(B-A)>.5,1); fixed=3-axisAlong;
lo=0; hi=1; ok=true;
if abs(D(fixed))>tol
    t=(A(fixed)-S(fixed))/D(fixed); q=S+t*D;
    ok=t>=-tol && t<=1+tol && q(axisAlong)>=min(A(axisAlong),B(axisAlong))-tol ...
        && q(axisAlong)<=max(A(axisAlong),B(axisAlong))+tol;
    lo=max(0,min(1,t)); hi=lo;
else
    if abs(S(fixed)-A(fixed))>tol, ok=false; return; end
    if abs(D(axisAlong))<=tol
        ok=S(axisAlong)>=min(A(axisAlong),B(axisAlong))-tol ...
            && S(axisAlong)<=max(A(axisAlong),B(axisAlong))+tol;
    else
        t=sort(([A(axisAlong) B(axisAlong)]-S(axisAlong))/D(axisAlong));
        lo=max(0,t(1)); hi=min(1,t(2)); ok=lo<=hi+tol;
    end
end
end

function Q=simplify(Q,tol)
j=2;
while j<=size(Q,1)
    if norm(Q(j,:)-Q(j-1,:))<=tol, Q(j,:)=[]; else, j=j+1; end
end
j=2;
while j<size(Q,1)
    a=Q(j,:)-Q(j-1,:); b=Q(j+1,:)-Q(j,:);
    if norm(cross(a,b))<=tol*(norm(a)+norm(b)) && dot(a,b)>=0
        Q(j,:)=[];
    else, j=j+1; end
end
end

function d=curveLength(Q)
d=sum(sqrt(sum(diff(Q,1,1).^2,2)));
end

function selftest()
cases={ [0.5 0 1],0,1; [0 0 1],.5,1; [0.5 0 0],1,1; ...
    [0.5 1 1],1,1; [0.5 1 0],2,2; [0 1 0],sqrt(3.25),2; ...
    [1 1 0],sqrt(3.25),2; [0.2 0 0.7],sqrt(.18),1; ...
    [1/3 1 0],sqrt(145)/6,4; [2/3 1 0],sqrt(145)/6,4};
for k=1:size(cases,1)
    p=cases{k,1}; [P,d,I]=cubeShortestPaths(p,false);
    assert(abs(d-cases{k,2})<1e-9,'Distance regression failed.');
    assert(numel(P)==cases{k,3},'Multiplicity regression failed.');
    assert(all(abs(I.pathLengths-d)<1e-9),'Polyline length failed.');
    fprintf('p=[%g %g %g]: length=%.12g, paths=%d\n',p,d,numel(P));
end

for x=[1/3-1e-6,1/3+1e-6]
    P=cubeShortestPaths([x 1 0],false); assert(numel(P)==2);
end
fprintf('All regression tests passed.\n');
end