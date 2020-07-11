function [Graph1,Graph2] = RipsCone(data1,data2,nv1,nv2,k)
nonZero1 = find(sum(data1,2));
nonZero2 = find(sum(data2,2));
newdata1=data1(nonZero1,:);
newdata2=data2(nonZero2,:);
nnv1=length(nonZero1);
nnv2=length(nonZero2);
newd1 = zeros(size(newdata1,1),2);
newd2 = zeros(size(newdata2,1),2);
r=10;
for i=1:size(newdata1,1)
    newd1(i,:) = plane2conic(newdata1(i,:),r);
end
for i=1:size(newdata2,1)
    newd2(i,:) = plane2conic(newdata2(i,:),r);
end
newdata1=newd1;
newdata2=newd2;
dis11 = zeros(size(newdata1,1));
dis22 = zeros(size(newdata2,1));
for i =1:size(newdata1,1)
    for j = 1:size(newdata1,1)
        dis11(i,j) = ConicDist(newdata1(i,:),newdata1(j,:));
    end
end
for i =1:size(newdata2,1)
    for j = 1:size(newdata2,1)
        dis22(i,j) = ConicDist(newdata2(i,:),newdata2(j,:));
    end
end
[~,id1] = sort(dis11,2);
[~,id2] = sort(dis22,2);
maxR=-1; 

G1=zeros(nnv1,nnv1);
G2=zeros(nnv2,nnv2);
for i=1:nnv1
    for j = i+1:nnv1
        ki = id1(i,2:k);
        kj = id1(j,2:k);
        if isempty(intersect(ki,kj)) == 0
            G1(i,j)=1;
            G1(j,i)=1;
        end
    end
end
for i=1:nnv2
    for j = i+1:nnv2
        ki = id2(i,2:k);
        kj = id2(j,2:k);
        if isempty(intersect(ki,kj)) == 0
            G2(i,j)=1;
            G2(j,i)=1;
        end
    end
end
if nnv1~=nv1
    newG1 = zeros(nv1,nnv1);
    k2=1;
    for k1=1:nv1
        if ~ismember(k1,nonZero1)
            newG1(k1,:)=0;
        else
            newG1(k1,:)=G1(k2,:);
            k2=k2+1;
        end
    end
    newG1 = [newG1 zeros(nv1, nv1-nnv1)];
    for k1=1:nv1
        if ~ismember(k1,nonZero1)
            newG1=[newG1(:,1:k1-1) circshift(newG1(:,k1:nv1),[0,1])];
        end
    end
    G1=newG1;
end
if nnv2~=nv2
    newG2 = zeros(nv2,nnv2);
    k2=1;
    for k1=1:nv2
        if ~ismember(k1,nonZero2)
            newG2(k1,:)=0;
        else
            newG2(k1,:)=G2(k2,:);
            k2=k2+1;
        end
    end
    newG2 = [newG2 zeros(nv2, nv2-nnv2)];
    for k1=1:nv2
        if ~ismember(k1,nonZero2)
            newG2=[newG2(:,1:k1-1) circshift(newG2(:,k1:nv2),[0,1])];
        end
    end
    G2=newG2;
end
R = corrcoef(G1,G2);
nz1 = find(sum(G1,2));
nz2 = find(sum(G2,2));
if R(1,2)>maxR && length(nz1)==length(nonZero1) && length(nz2)==length(nonZero2)
    maxR=R(1,2);
    Graph1=G1;
    Graph2=G2;
end
end
function c = plane2conic(coord,r)
x = coord(1);
y = coord(2);
% cartesian to polar on plane
R = sqrt(x^2 + y^2);
phi = atan2(y,x);
%polar on plane to spherical on sphere
Th = 2* atan(1/R);
theta = Th;
%spherical to 3d cartesian on sphere
xi = r * cos(phi*sin(theta));
si = r * sin(phi*sin(theta));
c = [xi,si];
end
function d = ConicDist(xi,xj)
d = sqrt((xi(1)-xj(1))^2 + (xi(2)-xj(2))^2);
end

function c = plane2ellipsoid(coord)
r=5;
x = coord(1);
y = coord(2);
R = sqrt(x^2 + y^2);
Th = atan2(y,x);
%polar on plane to spherical on sphere
phi = 2* atan(1/R);
theta = Th;
%spherical to 3d cartesian on sphere
x1 = r * cos(phi) * sin(theta);
y1 = r * sin(phi) * sin(theta);
z1 = r * cos(theta);
lon = atan2(y1,x1);
lat = acos(z1 / R);
c=[lat,lon];
end
function d = EllipsoidDist(lat1,lon1,lat2,lon2)
a = 6;
b=3;
f = (a-b)/a;
U1 = atan((1-f)*tan(lat1));
U2 = atan((1-f)*tan(lat2));
L = lon2-lon1;
lam=L;
sinsig = sqrt((cos(U2) *sin(lam))^2 + (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lam))^2);
cossig = sin(U1)*sin(U2) + cos(U1) * cos(U2) *cos(lam);
sig = atan(sinsig/cossig);
sinalpha = (cos(U1)*cos(U2)*sin(lam))/sinsig;
cos2alpha = 1- sinalpha.^2;
cos2sigm = cossig - ((2*sin(U1)*sin(U2))/cos2alpha);
C = (f/16)*cos2alpha*(4+f*(4-3*cos2alpha));
lam = L + (1-C)*f*sinalpha*(sig + C*sinsig*(cos2sigm + C*cossig * (-1+2*cos2sigm.^2)));
u2 = cos2alpha * ((a^2-b^2)/b^2);
A = 1 + (u2/16384)*(4096 + u2*(-768 + u2 * (320-175*u2)));
B = (u2/1024)*(256+u2*(-128+u2*(74-47*u2)));
dsig = B * sinalpha*(cos2sigm + (1/4)*B*(cossig*(-1+2*cos2sigm.^2) - (B/6) * cos2sigm * (-3+4*sinsig.^2)*(-3+4*cos2sigm.^2)));
d = b * A*(sig - dsig);
end

function c = plane2sphere(coord,r)
x = coord(1);
y = coord(2);
% cartesian to polar on plane
R = sqrt(x^2 + y^2);
Th = atan2(y,x);
%polar on plane to spherical on sphere
phi = 2* atan(1/R);
theta = Th;
%spherical to 3d cartesian on sphere
x1 = r * cos(phi) * sin(theta);
y1 = r * sin(phi) * sin(theta);
z1 = r * cos(theta);
c = [x1,y1,z1];
end
function d = sphericalDist(x,y,r)
chordal_distance = norm(x - y, 'fro');
d =  real(r*2*asin(.5*chordal_distance));
end
