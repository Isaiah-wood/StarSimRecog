close all
clear
format longG
%% 
StarLib = readmatrix('StarLib/hip_table.csv');

Fnsis = FuncStarImgSim();
Fnsr = FuncStarRecog();
Fnsm = FuncStarMatch();

[cameraConf, noiseConf] = Fnsis.InitConf();
%% 读取图像
[starImg, priorAtt, starList] = Fnsr.ReadStarImg(Fnsr.ImgDirPath,'ra63.3148_dec129.9164_roa170.455');
% [starImg, priorAtt, starList] = Fnsr.ReadStarImg(Fnsr.ImgDirPath);

%% 星图识别
binImg           = Fnsr.Binarization(starImg, 'timessigma', 3);
ObservedStarList = Fnsr.CenterExtraction(starImg, binImg, [], [], 3);

starVecListNum = 12;
FilteredStarList = Fnsm.StarListFilter(ObservedStarList, 12);

SortedStarCooList = FilteredStarList(:,2:-1:1);
SortedStarVecList =Fnsis.Coo2Vec(cameraConf, SortedStarCooList')';

starVecList = SortedStarVecList;
DistMatrix = zeros(starVecListNum);
for i = 1:starVecListNum
    Star1 = starVecList(i, :);
    for j = 1:starVecListNum
        if i == j
            DistMatrix(i, j) = 0;
            continue
        end
        Star2 = starVecList(j, :);
        DistMatrix(i, j) = dot(Star1, Star2)/(vecnorm(Star1)*vecnorm(Star2));
    end
end
[MaxDist, linear_Idx] = max(DistMatrix,[],"all");
[row, col] = ind2sub(size(DistMatrix), linear_Idx); 

poleStarList = [starVecList(row, :); starVecList(col, :)];

allStarList = SortedStarVecList;

Dist(:, :) = zeros(size(poleStarList, 1), size(allStarList, 1));
for i = 1:size(poleStarList, 1)
    poleStar = poleStarList(i, :);
    for j = 1:size(allStarList, 1)
        Star = allStarList(j, :);
        if poleStar == Star
            Dist(i, j) = 0;
            continue
        end
        Dist(i, j) = dot(poleStar, Star)/(vecnorm(poleStar)*vecnorm(Star));
    end

end
Distsorted = sort(Dist, 2, 'descend');

AngDist = Distsorted(1);

isBelowUpperBound = false(Fnsm.AngLibSize, size(allStarList, 1));
isAboveLowerBound = false(Fnsm.AngLibSize, size(allStarList, 1));
isBelowUpperBound(:,1) = any(AngDist > Fnsm.AngLibSubErr,2);
isAboveLowerBound(:,1) = any(AngDist < Fnsm.AngLibAddErr,2);
DistFilterMask = isBelowUpperBound & isAboveLowerBound;

MatchGroup = cell(size(poleStarList, 1), size(allStarList, 1));




MatchedSublib = Fnsm.AngLib(DistFilterMask(:,1),:);
stargroup = [MatchedSublib(:,2), MatchedSublib(:,3)];
MatchGroup{1,1} = stargroup;
MatchGroup{2,1} = stargroup;

for i = 1:size(poleStarList, 1)
    for j = 2:size(allStarList-1, 1)
        AngDist = Distsorted(i,j);
        isBelowUpperBound(:,j) = any(AngDist > Fnsm.AngLibSubErr,2);
        isAboveLowerBound(:,j) = any(AngDist < Fnsm.AngLibAddErr,2);
        DistFilterMask = isBelowUpperBound & isAboveLowerBound;
        MatchedSublib0 = Fnsm.AngLib(DistFilterMask(:,j),:);
        stargroup0 = [MatchedSublib0(:,2), MatchedSublib0(:,3)];
        MatchGroup{i,j} = intersect(MatchGroup{i,j-1},stargroup0);
    end
end

disp(['坐标为：', num2str(SortedStarCooList(row, :)), '的星点匹配到：'])
disp(MatchGroup{1,11})
disp(['坐标为：', num2str(SortedStarCooList(col, :)), '的星点匹配到：'])
disp(MatchGroup{2,11})
    

% 初始化逻辑数组，所有元素都为true  
keepColumns = true(size(MatchGroup, 2), 1);  
  
% 遍历MatchGroup的每一列  
for column = 1:size(MatchGroup, 2)  
    % 检查当前列是否包含{0×1 double}  
    if isempty(MatchGroup{1, column}) && size(MatchGroup{1, column}, 2) == 1  
        % 如果包含，将逻辑数组中对应位置的元素设置为false  
        keepColumns(column) = false;  
    end  
end  
  
% 使用逻辑数组来索引MatchGroup，并仅选择为true的列  
MatchGroupTrimmed = MatchGroup(:, keepColumns);  
  



% %% 计算姿态
% star1Id = MatchGroupTrimmed{1,end};
% star2Id = MatchGroupTrimmed{2,end};
% 
% 
% libVec0 = StarLib(StarLib(:,1) == star1Id,2:4);
% libVec1 = StarLib(StarLib(:,1) == star2Id,2:4);
% senVec0 = allStarList(row,:);
% senVec1 = allStarList(col,:);
% 
% 
% crossVec = cross(senVec0, senVec1);
% crossVec = crossVec / norm(crossVec);
% Ccm = [senVec0; crossVec; cross(senVec0, crossVec)];
% 
% crossVec = cross(libVec0, libVec1);
% crossVec = crossVec / norm(crossVec);
% Ccr = [libVec0; crossVec; cross(libVec0, crossVec)];
% 
% dcm = Ccm\Ccr;
% 
% att = Fnsis.Dcm2Att(dcm);
% disp(['计算得到像空间姿态：',num2str(rad2deg(att))])
% 
% priorDcm = Fnsis.Att2Dcm(deg2rad(priorAtt));
% calcAtt = Fnsis.Dcm2Att(priorDcm);
% Attr = deg2rad(priorAtt);

% disp([libVec0;libVec1;senVec0;senVec1])
% disp(Ccr/priorDcm)
% disp(Ccm)



% [starVecListbyDist, sortedStarListbyDist] = Fnsm.SortByAngDist(SortedStarVecList);
% starVecListbyDist(11,:) = SortedStarVecList(4,:);
% starGroup = [starVecListbyDist(1:end-1,:), sortedStarListbyDist, starVecListbyDist(2:end, :)];
% 
% satellite_num = size(starGroup, 1);
% MatchGroup               = cell(satellite_num,1);
% MatchGroupStarList       = cell(satellite_num,1);
% MatchGroupStarListSorted = cell(satellite_num,1);
% 
% isBelowUpperBound = false(Fnsm.AngLibSize, satellite_num);
% isAboveLowerBound = false(Fnsm.AngLibSize, satellite_num);
% for i = 1:size(starGroup, 1)
%     poleStar = starGroup(i, 1:3);
%     AngDist = starGroup(i, 4);
%     satellite = starGroup(i, 5:7);
%     isBelowUpperBound(:,i) = any(AngDist > Fnsm.AngLibSubErr,2);
%     isAboveLowerBound(:,i) = any(AngDist < Fnsm.AngLibAddErr,2);
%     DistFilterMask = isBelowUpperBound & isAboveLowerBound;
%     MatchGroup{i} = Fnsm.AngLib(DistFilterMask(:,i),:);
% 
%     MatchGroupStarList{i} = [MatchGroup{i}(:,2); MatchGroup{i}(:,3)];
%     MatchGroupStarListSorted{i} = unique(sort(MatchGroupStarList{i},'ascend'));
% 
% end


% att = [69.0283      132.9168       87.4259];
% dcm = Fnsis.Att2Dcm(deg2rad(att));
% 
% libVec0 = StarLib(1,2:4)';
% libVec1 = StarLib(2,2:4)';
% senVec0 = dcm * libVec0;
% senVec1 = dcm * libVec1;
% 
% 
% crossVec = cross(senVec0, senVec1);
% crossVec = crossVec / norm(crossVec);
% Ccm = [senVec0; crossVec; cross(senVec0, crossVec)];
% 
% crossVec = cross(libVec0, libVec1);
% crossVec = crossVec / norm(crossVec);
% Ccr = [libVec0; crossVec; cross(libVec0, crossVec)];
% 
% dcm_calc = Ccm\Ccr;
% 
% att_calc = Fnsis.Dcm2Att(dcm);
% 
% rad2deg(att_calc)


































% starGroup(:, 1) = FilteredStarList(:, 3);
% starGroup(:, 2:4) = starVecListbyDist;
% starGroup(:, 5) = FilteredStarList(:, 4);

% MatchedGroup1 = matchfunc(Fnsm, starGroup, 1);
% MatchedGroup2 = matchfunc(Fnsm, starGroup, 2, 1, MatchedGroup1);
% MatchedGroup3 = matchfunc(Fnsm, starGroup, 3, 2, MatchedGroup2);
% MatchedGroup4 = matchfunc(Fnsm, starGroup, 4, 3, MatchedGroup3);
% MatchedGroup5 = matchfunc(Fnsm, starGroup, 5, 4, MatchedGroup4);
% MatchedGroup6 = matchfunc(Fnsm, starGroup, 6, 5, MatchedGroup5);
% MatchedGroup7 = matchfunc(Fnsm, starGroup, 7, 6, MatchedGroup6);
% MatchedGroup8 = matchfunc(Fnsm, starGroup, 8, 7, MatchedGroup7);
% MatchedGroup9 = matchfunc(Fnsm, starGroup, 9, 8, MatchedGroup8);
% MatchedGroup10 = matchfunc(Fnsm, starGroup, 10, 9, MatchedGroup9);
% MatchedGroup11 = matchfunc(Fnsm, starGroup, 11, 10, MatchedGroup10);










%% 函数
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

        cosTheta = dot(poleStar, satVector)/(vecnorm(poleStar)*vecnorm(satVector));
        AngDist(i) = cosTheta;
    end
end