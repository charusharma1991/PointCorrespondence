function [H, neighborhood, neighborCoeff] = Ellipsoidal_Weights( pts, E )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the reconstruction weights for each point in a
% point-set.
% 
% It should be called by Convex_Matching_Affine.m 
%
% INPUT
% pts           An N x 2 matrix which represents a N point-set
% E             An N x N 0-1 edge matrix. The ith row records the ith 
%               point's neighbors
% 
% OUTPUT
% H             The N x N reconstruction matrix. H = I - W. 
%               Ideally, H * pts should be an N x 2 zero matrix.
% neighborhood  An N x 1 cell vector. The ith element records the ith point's 
%               neighbors' indices.
% neighborCoeff An N x 1 cell vector. The ith element records the ith point's 
%               neighbors' reconstruction weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nm, d] = size( pts );

neighborCoeff = cell( Nm, 1 );
neighborhood = cell( Nm, 1 );

for i = 1:Nm
    neighborhood{i} = find( E(i,:) ~= 0 );
    
    if length( neighborhood{i} ) < 1
        error('E matrix illegal!\n The %dth point has at most 0 neighbors!\n', i);
    end
end

% Calculate the convex combination
H = sparse( Nm, Nm );
for i = 1:Nm
    
     % current point
     m = pts(i,:);
     m = plane2ellipsoid(m);
     % all neighbors
     n = pts(neighborhood{i}, : );
     
     % number of three neighbors
     nNum = length( neighborhood{i} );
     nn=[];
     dv=[];
     for j = 1:nNum
         p = plane2ellipsoid(n(j,:));
         nn = [nn; p];
         %disp(p)
         d = EllipsoidDist(m(1),m(2),p(1),p(2));
         dv = [dv,d];
     end
     NormRows = sqrt(sum(dv.*dv,2));
     Ynorm = bsxfun(@rdivide,abs(dv),NormRows);
     dv = Ynorm;

     neighborCoeff{i} = dv;
     H( i, neighborhood{i} ) = dv;
     H( i, i ) = 1;
end
H = sparse(H);
end
function c = plane2ellipsoid(coord)
r=10;
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
c = [lat,lon];
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
