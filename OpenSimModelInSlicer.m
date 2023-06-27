%-------------------------------------------------------------------------%
% Copyright (c) 2022 Koller W.                                            %
%    Author:   Willi Koller,  2022                                        %
%    email:    willi.koller@univie.ac.at                                  %
% ----------------------------------------------------------------------- %

% Method_Sangeux_2015 (c) by Willi Koller, University of Vienna
%
% Method_Sangeux_2015 is licensed under a
% Creative Commons Attribution-NonCommercial 4.0 International License.
%
% You should have received a copy of the license along with this
% work. If not, see <http://creativecommons.org/licenses/by-nc/4.0/>.

% this script uses GIBBON- and MSK-STAPLE toolboxes
% https://github.com/modenaxe/msk-STAPLE
% https://github.com/gibbonCode/

clear;
addpath(genpath('./lib'));
% STLfileName = [];

% Kira Right Tibia
fname = 'C:\Users\Willi\ucloud\BaseAngles\HamnerGeometryInSlicer\R_Fiducials_Hamner.mrk.json';
side = 'R';
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
data = jsondecode(str);
% Tibia torsion in 3D = 27.7°
% Tibial Torsion in 2D = 26.83°
% Neck-Shaft Angle = 119.3°
% Anteversion in plane perpendicular to shaft = 18.3°
% Anteversion in transverse plane (OpenSim coordinate system) = 17.4°

% % Kira Left Tibia
% fname = 'C:\Users\Willi\ucloud\BaseAngles\HamnerGeometryInSlicer\L_Fiducials_Hamner.mrk.json';
% side = 'L';
% fid = fopen(fname);
% raw = fread(fid,inf);
% str = char(raw');
% fclose(fid);
% data = jsondecode(str);
% % Tibia torsion in 3D = 23.8°
% % Tibial Torsion in 2D = 23.82°
% % Neck-Shaft Angle = 119°
% % Anteversion in plane perpendicular to shaft = 19.2°
% % Anteversion in transverse plane (OpenSim coordinate system) = 18.3°


for i = 1 : length(data.markups.controlPoints)
    switch data.markups.controlPoints(i).label
        case 'MT'
            mtI = i;
        case 'LT'
            ltI = i;
        case 'MMAL'
            mmalI = i;
        case 'LMAL'
            lmalI = i;        
    end
end

tibiaPlateauAxis = data.markups.controlPoints(mtI).position - data.markups.controlPoints(ltI).position;
malleoliAxis = data.markups.controlPoints(mmalI).position - data.markups.controlPoints(lmalI).position;

TibiaTorsion = atan2(norm(cross(tibiaPlateauAxis, malleoliAxis)), dot(tibiaPlateauAxis, malleoliAxis));
TibiaTorsion_deg = TibiaTorsion*180/pi;

if TibiaTorsion_deg > 90
    TibiaTorsion_deg = 180-TibiaTorsion_deg;
end
disp(['Tibia torsion in 3D = ' num2str(TibiaTorsion_deg, 3) '°']);

% kira code
tibialAxis = tibiaPlateauAxis;
ankleAxis = malleoliAxis;

% TT in 2D
u = [1;0;0];
v = [0;0;1];
% % u = [0;0;1];
% % v = [1;0;0];
tibialAxis = dot(tibialAxis,u)*u + dot(tibialAxis,v)*v; % orthogonal projection onto the x-y-plane
ankleAxis = dot(ankleAxis,u)*u + dot(ankleAxis,v)*v;
TT = atan2d(vecnorm(cross(tibialAxis,ankleAxis)),dot(tibialAxis,ankleAxis));

if TT > 90
    TT = 180 - TT;
end
disp(['Tibial Torsion in 2D = ' num2str(TT, 4) '°']);

% FEMUR
for i = 1 : length(data.markups.controlPoints)
    switch data.markups.controlPoints(i).label
        case 'HJC'
            hjcI = i;
        case 'GT'
            gtI = i;
        case 'PS'
            psI = i;
        case 'DS'
            dsI = i;
        case 'LC'
            lcI = i;
        case 'MC'
            mcI = i;
    end
end

neckAxis = data.markups.controlPoints(hjcI).position - data.markups.controlPoints(gtI).position;
shaftAxis = data.markups.controlPoints(psI).position - data.markups.controlPoints(dsI).position;
kneeAxis = data.markups.controlPoints(lcI).position - data.markups.controlPoints(mcI).position;


% NSA in 3D
NSA = atan2d(vecnorm(cross(neckAxis*(-1),shaftAxis)),dot(neckAxis*(-1),shaftAxis));
disp(['Neck-Shaft Angle = ' num2str(NSA, 4) '°']);

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0; ...
              norm(cross(A,B)) dot(A,B)  0; ...
              0              0           1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

UU = @(Fi,G) Fi*G*inv(Fi);
a=shaftAxis; 
a = a /norm(a);
b=[0 0 1]';
U = UU(FFi(a,b), GG(a,b));
norm(U);
norm(b-U*a);
U;

neckAxisRotated = U * neckAxis;
kneeAxisRotated = U * kneeAxis;
shaftAxisRotated = U * shaftAxis;

u = neckAxisRotated([1 2]);
v = kneeAxisRotated([1 2]) *(-1);
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
AVA_inPlanePerpendicularToShaft = real(acosd(CosTheta));
disp(['Anteversion in plane perpendicular to shaft = ' num2str(AVA_inPlanePerpendicularToShaft, 3) '°']);

figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.5]);
tiledlayout(1, 2);
nexttile;
hold on;
rotate3d on;
daspect([1 1 1]);
quiver3(data.markups.controlPoints(gtI).position(1) ,data.markups.controlPoints(gtI).position(2) ,data.markups.controlPoints(gtI).position(3) , neckAxis(1), neckAxis(2), neckAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(dsI).position(1) ,data.markups.controlPoints(dsI).position(2) ,data.markups.controlPoints(dsI).position(3) , shaftAxis(1), shaftAxis(2), shaftAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(mcI).position(1) ,data.markups.controlPoints(mcI).position(2) ,data.markups.controlPoints(mcI).position(3) , kneeAxis(1), kneeAxis(2), kneeAxis(3), 'LineWidth', 2);
leg = legend({'Neck Axis', 'Shaft Axis', 'Knee Axis'});
leg.Layout.Tile = "north";
title('View on plane which is perpendicular to shaft axis')
view(shaftAxis);

nexttile;
hold on;
rotate3d on;
daspect([1 1 1]);
quiver3(data.markups.controlPoints(gtI).position(1) ,data.markups.controlPoints(gtI).position(2) ,data.markups.controlPoints(gtI).position(3) , neckAxis(1), neckAxis(2), neckAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(dsI).position(1) ,data.markups.controlPoints(dsI).position(2) ,data.markups.controlPoints(dsI).position(3) , shaftAxis(1), shaftAxis(2), shaftAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(mcI).position(1) ,data.markups.controlPoints(mcI).position(2) ,data.markups.controlPoints(mcI).position(3) , kneeAxis(1), kneeAxis(2), kneeAxis(3), 'LineWidth', 2);
title('View in transverse plane (if in OpenSim coordinate system)')
view(0, 0);

u = neckAxis([1 3]);
v = kneeAxis([1 3]) *(-1);
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
AVA_inTransverse = real(acosd(CosTheta));
disp(['Anteversion in transverse plane (OpenSim coordinate system) = ' num2str(AVA_inTransverse, 3) '°']);