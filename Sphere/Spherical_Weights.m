function [H, neighborhood, neighborCoeff] = Spherical_Weights( pts, E )

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
% H = I - W
H = sparse( Nm, Nm );
r=5;
for i = 1:Nm
     % current point
     m = pts(i,:);
     m = plane2sphere(m,r);
     % all neighbors
     n = pts(neighborhood{i}, : );
     
     % number of three neighbors
     nNum = length( neighborhood{i} );
     nn=[];
     dv=[];
     for j = 1:nNum
         p = plane2sphere(n(j,:),r);
         nn = [nn; p];
         d = sphericalDist(m,p,r);
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
