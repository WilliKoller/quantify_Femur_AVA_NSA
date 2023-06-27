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
STLfileName = [];

% Hans Right Tibia
side = 'R';


% Hans Left Tibia
side = 'L';


% Willi Right Tibia
side = 'R';
data.markups.controlPoints(1).label = 'MT';
data.markups.controlPoints(1).position = [-0.031363174799423001 -0.035767237677210233 -0.010505321939548623]';
data.markups.controlPoints(2).label = 'LT';
data.markups.controlPoints(2).position = [-0.028637666866231447 -0.032512271155183896 0.024452631058088769]';
data.markups.controlPoints(3).label = 'MMAL';
data.markups.controlPoints(3).position = [0.0048342479309876228 -0.37966037290303201 -0.026218792674424479]';
data.markups.controlPoints(4).label = 'LMAL';
data.markups.controlPoints(4).position = [-0.016858826944702095 -0.38602392366675503 0.036037096620815105]';
% Tibia torsion in 3D = 26°
% Tibial Torsion in 2D = 33.71°

% Willi Left Tibia
side = 'L';
data.markups.controlPoints(1).label = 'MT';
data.markups.controlPoints(1).position = [-0.031056442873601788 -0.035767237677210212 0.0085228399613840877]';
data.markups.controlPoints(2).label = 'LT';
data.markups.controlPoints(2).position = [-0.029643970452373002 -0.032512271155183875 -0.020019428887615598]';
data.markups.controlPoints(3).label = 'MMAL';
data.markups.controlPoints(3).position = [0.0022075017976223251 -0.37829281200248122 0.028332522929387066]';
data.markups.controlPoints(4).label = 'LMAL';
data.markups.controlPoints(4).position = [-0.017344969569836013 -0.38460564098538935 -0.038043035005905473]';
% Tibia torsion in 3D = 22.5°
% Tibial Torsion in 2D = 48.65°


% % Willi Right Tibia
% fname = 'right_tibia_willi.mrk.json';
% side = 'R';
% fid = fopen(fname);
% raw = fread(fid,inf);
% str = char(raw');
% fclose(fid);
% data = jsondecode(str);
% % Tibia torsion in 3D = 25.2°
% % Tibial Torsion in 2D = 18.71°


% % Willi Left Tibia
% fname = 'left_tibia_willi.mrk.json';
% side = 'L';
% fid = fopen(fname);
% raw = fread(fid,inf);
% str = char(raw');
% fclose(fid);
% data = jsondecode(str);
% 
% % Tibia torsion in 3D = 26.4°
% % Tibial Torsion in 2D = 8.46°
 

% % Bas Right Tibia
% fname = 'points_r_Bas.mrk.json';
% side = 'R';
% fid = fopen(fname);
% raw = fread(fid,inf);
% str = char(raw');
% fclose(fid);
% data = jsondecode(str);
% % Tibia torsion in 3D = 19.5°
% % Tibial Torsion in 2D = 47.92°
% 
% 
% % Bas Left Tibia
% fname = 'points_l_Bas.mrk.json';
% side = 'L';
% fid = fopen(fname);
% raw = fread(fid,inf);
% str = char(raw');
% fclose(fid);
% data = jsondecode(str);
% % Tibia torsion in 3D = 11.4°
% % Tibial Torsion in 2D = 21.83°

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
v = [0;1;0];
% % u = [0;0;1];
% % v = [1;0;0];
tibialAxis = dot(tibialAxis,u)*u + dot(tibialAxis,v)*v; % orthogonal projection onto the x-y-plane
ankleAxis = dot(ankleAxis,u)*u + dot(ankleAxis,v)*v;
TT = atan2d(vecnorm(cross(tibialAxis,ankleAxis)),dot(tibialAxis,ankleAxis));

if TT > 90
    TT = 180 - TT;
end
disp(['Tibial Torsion in 2D = ' num2str(TT, 4) '°']);
