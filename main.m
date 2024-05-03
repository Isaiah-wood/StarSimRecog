close all
clear

%% 引用
sif = StarSensorInitFunc();
isf = StarImgSimFunc();
srf = StarRecogFunc();


%% 模拟生成星图
[sensorConf, ~, imgBackgdConf, noiseConf] = isf.ExampleConf();
% 通过改sensorConf、imgBackgdConf、noiseConf中的值来修改星图模拟参数
imgPointConf = isf.GetStar(sif, sensorConf, sif.StarLib, isf.StarLibInVision);
starImg = isf.ImgGen(sensorConf, imgPointConf, imgBackgdConf, noiseConf);

% imshow(starImg);
% isf.SaveImg(sensorConf,starImg);

%% 星图识别
binImg = srf.Binarization(starImg, 'TimesSigma', 2.2);
starList = srf.CenterExtraction(starImg, binImg, [], [], 3);
markImg = srf.MarkPosition(starImg,starList);

imshow(markImg);



%%
[sensorConf.ra, sensorConf.dec, sensorConf.roa]
att