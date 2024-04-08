close all
clear

%% 引用
sif = StarSensorInitFunc();
isf = StarImgGenFunc();


%% 模拟生成星图
[sensorConf, ~, imgBackgdConf, noiseConf] = isf.ExampleConf();
% 通过改sensorConf、imgBackgdConf、noiseConf中的值来修改星图模拟参数
imgPointConf = isf.GetStar(sif, sensorConf, sif.starLib60, isf.StarLibParserSAO60SSS);
starImg = isf.ImgGen(sensorConf, imgPointConf, imgBackgdConf, noiseConf);
imshow(starImg);