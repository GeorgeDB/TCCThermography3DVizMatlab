% =========================================================================
% Universidade Federal de Uberl�ndia - UFU
% C�digo TCC - Fus�o multimodal para visualiza��o 3D de termografia
% infravemelha
% Autor: George Dechichi Barbar - 11721EMC020
% =========================================================================

%%
clear; clc; close all;
%%

% Nomes dos arquivos a serem importados
%%
cameraCalibrationFileName = 'calibSession.mat';
modelFileName = '3dModel.stl';
thermographyFileName = 'thermographyResults.mat';

load(cameraCalibrationFileName);
load(thermographyFileName);
%%

% Gera��o da malha no matlab
%%
numberOfPDEs = 1;
pdem = createpde(numberOfPDEs);
gm = importGeometry(pdem,modelFileName);
msh = generateMesh(pdem,'Hmax',1,'Hmin',0.01,'Hgrad',2);
malha = msh.Nodes';
%%

% Matrizes de rota��o e transla��o
%%
imagePoints = [];
worldPoints = [];
focalLength = mean(calibrationSession.CameraParameters.FocalLength);
intrinsicMatrix = calibrationSession.CameraParameters.IntrinsicMatrix;
center = [intrinsicMatrix(3,1) intrinsicMatrix(3,2)];
[rot,trans] = modernPosit(imagePoints, worldPoints, focalLength, center);
%%

% Pontos na refer�ncia da c�mera
%%
meshPoints = msh.Nodes;
meshPoints(4,:) = 1;
world2CameraMatrix = [rot trans];
world2CameraMatrix(4,4) = 1;
world2Camera = world2CameraMatrix * meshPoints;
%%

% Proje��o na percep��o da c�mera
%%
cameraPerspectiveMatrix = instrinsicMatrix';
cameraPerspectiveMatrix(:,4) = 0;
worldProjection = cameraPerspectiveMatrix * world2Camera;
%%

% Rela��o entre pixels e pontos
%%
x_linha = worldProjection(1,:);
y_linha = worldProjection(2,:);
z_linha = worldProjection(3,:);
numUVPixels = length(x_linha);
pixelsMatrix = zeros(numUVPixels,2);
temperatureMesh = zeros(numUVPixels,4);

frame = 160;
for i = 1:1:numUVPixels
    u = round(x_linha(i)/z_linha(i));
    v = round(y_linha(i)/z_linha(i));
    T_ref = Raw(u, v, frame);
    pixelsMatrix(i,:) = [u v];
    temperatureMesh(i,:) = [malha(i,:) T_ref];
end
%%

% Plot para conferir semelhan�a de pontos
%%
sizeRawImage = size(Raw);
imagesc(Raw(1:sizeRawImage(1,1),1:sizeRawImage(1,2),frame));
hold on;
plot(pixelsMatrix(:,1),pixelsMatrix(:,2),'.k')
%%



