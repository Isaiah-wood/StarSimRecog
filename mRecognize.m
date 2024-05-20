close all
clear

%% 
StarLib = readmatrix('StarLib/hip_table.csv');

Fnsis = FuncStarImgSim();
Fnsr = FuncStarRecog();
Fnsm = FuncStarMatch();

[cameraConf, noiseConf] = Fnsis.InitConf();
%% 读取图像
starImg = Fnsr.ReadStarImg(Fnsr.ImgDirPath,'ra112.5788_dec29.0673_roa64.3558.png');

%% 星图识别
binImg = Fnsr.Binarization(starImg, 'TimesSigma', 3);
starList = Fnsr.CenterExtraction(starImg, binImg, [], [], 3);
FilteredStarList = Fnsm.StarListFilter(starList, 12); 





% FilteredStarListSt.bri = FilteredStarList(:, 1);
% FilteredStarListSt.row = FilteredStarList(:, 2);
% FilteredStarListSt.col = FilteredStarList(:, 3);
% FilteredStarListSt.size = FilteredStarList(:, 4);
% markImg = Fnsr.MarkPosition(starImg,FilteredStarListSt);
% imshow(markImg);
% imwrite(markImg,'topbristar.png');
% starGroup = sortrows(FilteredStarList,1,"descend");
% [distMat,~] = Fnsm.AngMatch(Fnsm.AngLib, starGroup, 12);



SortedStarList = sortrows(FilteredStarList,1,"descend");

SortedStarCooList = SortedStarList(:,2:3);
SortedStarVecList =Fnsis.Coo2Vec(cameraConf, SortedStarCooList')';

starGroup(:, 1) = FilteredStarList(:, 1);
starGroup(:, 2:4) = SortedStarVecList;
starGroup(:, 5) = FilteredStarList(:, 4);

MatchedGroup1 = matchfunc(Fnsm, starGroup, 1);
MatchedGroup2 = matchfunc(Fnsm, starGroup, 2, 1, MatchedGroup1);
MatchedGroup3 = matchfunc(Fnsm, starGroup, 3, 2, MatchedGroup2);
MatchedGroup4 = matchfunc(Fnsm, starGroup, 4, 3, MatchedGroup3);
MatchedGroup5 = matchfunc(Fnsm, starGroup, 5, 4, MatchedGroup4);
MatchedGroup6 = matchfunc(Fnsm, starGroup, 6, 5, MatchedGroup5);
MatchedGroup7 = matchfunc(Fnsm, starGroup, 7, 6, MatchedGroup6);
MatchedGroup8 = matchfunc(Fnsm, starGroup, 8, 7, MatchedGroup7);
MatchedGroup9 = matchfunc(Fnsm, starGroup, 9, 8, MatchedGroup8);
MatchedGroup10 = matchfunc(Fnsm, starGroup, 10, 9, MatchedGroup9);
MatchedGroup11 = matchfunc(Fnsm, starGroup, 11, 10, MatchedGroup10);
MatchedGroup1 = matchfunc(Fnsm, starGroup, 1, 11, MatchedGroup11);
MatchedGroup2 = matchfunc(Fnsm, starGroup, 2, 1, MatchedGroup1);
MatchedGroup3 = matchfunc(Fnsm, starGroup, 3, 2, MatchedGroup2);
MatchedGroup4 = matchfunc(Fnsm, starGroup, 4, 3, MatchedGroup3);
MatchedGroup5 = matchfunc(Fnsm, starGroup, 5, 4, MatchedGroup4);
MatchedGroup6 = matchfunc(Fnsm, starGroup, 6, 5, MatchedGroup5);
MatchedGroup7 = matchfunc(Fnsm, starGroup, 7, 6, MatchedGroup6);
MatchedGroup8 = matchfunc(Fnsm, starGroup, 8, 7, MatchedGroup7);
MatchedGroup9 = matchfunc(Fnsm, starGroup, 9, 8, MatchedGroup8);
MatchedGroup10 = matchfunc(Fnsm, starGroup, 10, 9, MatchedGroup9);
MatchedGroup11 = matchfunc(Fnsm, starGroup, 11, 10, MatchedGroup10);
MatchedGroup = [MatchedGroup1, MatchedGroup2, MatchedGroup3, MatchedGroup4, MatchedGroup5, ...
    MatchedGroup6, MatchedGroup7, MatchedGroup8, MatchedGroup9, MatchedGroup10, MatchedGroup11]';
% disp(MatchedGroup)

% MatchedGroup1 = matchfunc(Fnsm, starGroup, 5);
% MatchedGroup2 = matchfunc(Fnsm, starGroup, 4, 5, MatchedGroup1);
% MatchedGroup3 = matchfunc(Fnsm, starGroup, 3, 4, MatchedGroup2);
% MatchedGroup4 = matchfunc(Fnsm, starGroup, 2, 3, MatchedGroup3);
% MatchedGroup5 = matchfunc(Fnsm, starGroup, 1, 2, MatchedGroup4);

star1Id = find(StarLib(:,1) == MatchedGroup1);
star2Id = find(StarLib(:,1) == MatchedGroup2);


libVec0 = StarLib(star1Id,2:4);
libVec1 = StarLib(star2Id,2:4);
senVec0 = SortedStarVecList(1,:);
senVec1 = SortedStarVecList(2,:);
senVec0 = flip(senVec0);
senVec1 = flip(senVec1);

dcm = Fnsm.DualVecAttiCalc(libVec0, libVec1, senVec0, senVec1);
att = Fnsis.Dcm2Att(dcm);
disp(att)


att_ = [3.8491, 2.1456, 4.8080];
dcm_ = Fnsis.Att2Dcm(att_);

% crosssenVec_ = cross(senVec0, senVec1);
% crosssenVec_ = crosssenVec_ / norm(crosssenVec_);
% Ccm_ = [senVec0; crosssenVec_; cross(senVec0, crosssenVec_)];

crosslibVec_ = cross(libVec0, libVec1);
crosslibVec_ = crosslibVec_ / norm(crosslibVec_);
Ccr_ = [libVec0; crosslibVec_; cross(libVec0, crosslibVec_)];


Ccm_ = Ccr_/dcm;





%%
function MatchedGroup = matchfunc(Fnsm, starGroup, poleStarSn, prepoleStarSn, preMatchedGroup)
    poleStar = starGroup(poleStarSn, 2:4);
    satelliteSn = 1:size(starGroup, 1);
    satelliteSn(poleStarSn) = [];
    satellite = starGroup(satelliteSn, 2:4);


    if nargin <= 3
        preMatchedFlag = 0;
    else
        preMatchedFlag = 1;
        if prepoleStarSn < poleStarSn
            preMatchedStarSn = prepoleStarSn;
        elseif prepoleStarSn > poleStarSn
            preMatchedStarSn = prepoleStarSn - 1;
        else
            error ('PreMatchedStarSn and prepoleStarSn should not be the same')
        end
    end
    satellite_num = size(satellite, 1);
    MatchGroup = cell(satellite_num,1);
    MatchGroupStarList = cell(satellite_num,1);
    MatchGroupStarListSorted = cell(satellite_num,1);

    AngDist = GetAngDist(poleStar, satellite);
    isBelowUpperBound = false(Fnsm.AngLibSize, satellite_num);
    isAboveLowerBound = false(Fnsm.AngLibSize, satellite_num);
    for i = 1:satellite_num
        isBelowUpperBound(:,i) = any(AngDist(i) > Fnsm.AngLibSubErr,2);
        isAboveLowerBound(:,i) = any(AngDist(i) < Fnsm.AngLibAddErr,2);
        DistFilterMask = isBelowUpperBound & isAboveLowerBound;

        MatchGroup{i} = Fnsm.AngLib(DistFilterMask(:,i),:);

        if preMatchedFlag == 0
            MatchGroupStarList{i} = [MatchGroup{i}(:,2); MatchGroup{i}(:,3)];
            MatchGroupStarListSorted{i} = unique(sort(MatchGroupStarList{i},'ascend'));
        elseif preMatchedFlag == 1
            if i == preMatchedStarSn
                for j = size(MatchGroup{i},1):-1:1
                    if ismember(MatchGroup{i}(j,2), preMatchedGroup) 
                        MatchGroup{i}(j,2) = MatchGroup{i}(j,3);
                    elseif ismember(MatchGroup{i}(j,3), preMatchedGroup)
                        MatchGroup{i}(j,3) = MatchGroup{i}(j,2);
                    else
                        MatchGroup{i}(j,:) = [];
                    end
                end
            end
            MatchGroupStarList{i} = [MatchGroup{i}(:,2); MatchGroup{i}(:,3)];
            MatchGroupStarListSorted{i} = unique(sort(MatchGroupStarList{i},'ascend'));
        end



        if i == 1 
            common_elements = MatchGroupStarListSorted{1};
        else 
            common_elements = intersect(common_elements, MatchGroupStarListSorted{i});
        end
        % disp(['序号：',num2str(i)])
        % disp(size(common_elements))
    end
    % disp(common_elements)

    MatchedGroup = common_elements;
end


function AngDist = GetAngDist(poleStar, satellite)
    % 检查输入的维度
    assert(size(poleStar, 2) == 3, 'poleStar must be a 1x3 vector');
    assert(all(size(satellite, 2) == 3), 'All rows of satellite must be 3-dimensional vectors');

    % 规范化输入的向量
    poleStarNorm = norm(poleStar);
    poleStar = poleStar / poleStarNorm;
    
    % 初始化结果数组
    AngDist = zeros(size(satellite, 1), 1);

    % 计算夹角并转换为度
    for i = 1:size(satellite, 1)
        satVector = satellite(i, :);  % 获取卫星向量
        satVectorNorm = norm(satVector);
        satVector = satVector / satVectorNorm;

        cosTheta = dot(poleStar, satVector);
        AngDist(i) = acosd(cosTheta);
    end
end