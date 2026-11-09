function polygon = validatePsiPolygon(polygon)
%VALIDATEPSIPOLYGON Validate a simple geographic ring; unwrap the date line.
if ~isnumeric(polygon) || size(polygon,2)~=2 || size(polygon,1)<3 || ...
        any(~isfinite(polygon),'all') || any(abs(polygon(:,1))>180) || ...
        any(abs(polygon(:,2))>=90)
    error('PHASE:InvalidPsiArea','Draw a polygon with at least three valid longitude/latitude vertices.');
end
polygon=double(polygon);
if ~isequal(polygon(1,:),polygon(end,:)),polygon(end+1,:)=polygon(1,:);end
p=polygon;p(:,1)=rad2deg(unwrap(deg2rad(p(:,1))));
if max(p(:,1))-min(p(:,1))>=180 || polyarea(p(:,1),p(:,2))<1e-12 || ...
        size(unique(p(1:end-1,:),'rows'),1)~=size(p,1)-1
    error('PHASE:InvalidPsiArea','The PSI polygon must have distinct vertices and a nonzero area smaller than a hemisphere.');
end
n=size(p,1)-1;
for i=1:n
    for j=i+1:n
        if j==i+1 || (i==1 && j==n),continue;end
        a=p(i,:);b=p(i+1,:);c=p(j,:);d=p(j+1,:);
        if intersectSegments(a,b,c,d)
            error('PHASE:InvalidPsiArea','The PSI polygon crosses or touches itself. Redraw its boundary.');
        end
    end
end
end

function yes=intersectSegments(a,b,c,d)
cross2=@(u,v)u(1)*v(2)-u(2)*v(1);
o=[cross2(b-a,c-a),cross2(b-a,d-a),cross2(d-c,a-c),cross2(d-c,b-c)];
tol=1e-12;
yes=o(1)*o(2)<-tol && o(3)*o(4)<-tol;
p={c,d,a,b};lo={min(a,b),min(a,b),min(c,d),min(c,d)};
hi={max(a,b),max(a,b),max(c,d),max(c,d)};
for k=1:4
    yes=yes || (abs(o(k))<=tol && all(p{k}>=lo{k}-tol) && all(p{k}<=hi{k}+tol));
end
end
