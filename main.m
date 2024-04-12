close all
clear

%% 引用
sif = StarSensorInitFunc();
isf = StarImgSimFunc();


%% 模拟生成星图
[sensorConf, ~, imgBackgdConf, noiseConf] = isf.ExampleConf();
% 通过改sensorConf、imgBackgdConf、noiseConf中的值来修改星图模拟参数
imgPointConf = isf.GetStar(sif, sensorConf, sif.StarLib, isf.StarLibInVision);
starImg = isf.ImgGen(sensorConf, imgPointConf, imgBackgdConf, noiseConf);
imshow(starImg);