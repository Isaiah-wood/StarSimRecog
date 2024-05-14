close all
clear


StarLib = readmatrix('StarLib/hip_table.csv');
Fnsis = FuncStarImgSim();

[cameraConf, noiseConf] = Fnsis.InitConf();
VisibleStarList = Fnsis.StarLibInVision(StarLib, cameraConf);
imgPointConf    = Fnsis.TakePhoto(cameraConf, VisibleStarList);
starImg         = Fnsis.PrintPhoto(cameraConf, imgPointConf, noiseConf);

imshow(starImg);
% Fnsis.SaveImgWithDir(cameraConf,starImg);

% VisibleStarListSorted = sortrows(VisibleStarList,5);
% writematrix(VisibleStarListSorted,'VisibleStarListSorted.csv')