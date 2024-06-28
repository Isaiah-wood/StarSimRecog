function func = FuncStarMatch
    func.WorkspacePath = strrep(fileparts(mfilename('fullpath')), '\', '/');
    % func.AngLibPath = [func.WorkspacePath, '/StarLib/angle_database.csv'];
    func.AngLibPath = [func.WorkspacePath, '/StarLib/AngLib_mag8_sorted.csv'];
    func.AngLib = readmatrix(func.AngLibPath);
    func.AngLibSize = size(func.AngLib, 1);

    func.matchconfig.starNumMax = 12;
    func.matchconfig.errordist = 1e-07;
    func.matchconfig.errordbri = 0.05;

    func.AngLibSubErr = func.AngLib(:, 1) - func.matchconfig.errordist;
    func.AngLibAddErr = func.AngLib(:, 1) + func.matchconfig.errordist;


    func.AngMatch = @AngMatch;
    func.StarListFilter = @StarListFilter;
    func.SortByAngDist = @SortByAngDist;
    func.DualVecAttiCalc = @DualVecAttiCalc;
end

function FilteredStarList = StarListFilter(starList, starNumMax)
    bricol = 4;
    sortedBriList = sortrows(starList, bricol, 'descend');
    % briShel = sortedBriList(starNumMax+1, bricol);
    FilteredStarList(:, 1:4) = sortedBriList(1:starNumMax,:);
end


function [starVecListbyDist, sortedStarListbyDist] = SortByAngDist(starVecList)
    starVecListNum = size(starVecList, 1);
    DistMatrix = zeros(starVecListNum);

    starVecListbyDist = zeros(starVecListNum, 3);
    sortedStarListbyDist = zeros(starVecListNum - 1, 1);
    % starVecListCopy = starVecList;

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
    % disp(DistMatrix)
    [MaxDist, linear_Idx] = max(DistMatrix,[],"all");
    % disp(MaxDist)
    [row, col] = ind2sub(size(DistMatrix), linear_Idx); 
    % disp([row, col])
    sortedStarListbyDist(1) = MaxDist(1);
    starVecListbyDist(1,:) = starVecList(row, :);



    DistMatrix(row, :) = 0;
    % disp(DistMatrix)
    poleStarid = col;
    for i = 2:starVecListNum
        % poleStar = starVecList(poleStarid,:);
        [MaxDist, linear_Idx] = max(DistMatrix(:, poleStarid));
        % disp(linear_Idx)
        % [row, col] = ind2sub(size(DistMatrix), linear_Idx); 
        sortedStarListbyDist(i) = MaxDist;
        starVecListbyDist(i,:) = starVecList(poleStarid, :);
        poleStarid = linear_Idx;
        DistMatrix(linear_Idx, :) = 0;
    disp(DistMatrix)
    end
    sortedStarListbyDist(12) = [];
    

    % CloseStarGroup = [starVecList(row, :),starVecList(col, :)];
    % Distlist = [MaxDist(row, 1:col-1),0,MaxDist(row, col+1:end);...
    %             MaxDist(1:row-1, col),0,MaxDist(row+1:end, col)];
    % [MaxDist, linear_Idx] = max(Distlist);
    % [row, col] = ind2sub(size(Distlist), linear_Idx); 
    % sortedStarListbyDist(2) = MaxDist;
    % if row == 1 
    %     starVecListbyDist(1:2, :) = CloseStarGroup(2:1);
    % elseif row == 2 
    %     starVecListbyDist(1:2, :) = CloseStarGroup(1:2);
    % end
    % starVecListbyDist(3, :) = starVecList(col, :);
    % Distlist = MaxDist(col, :);
    





end










%% 匹配组算法
% 1. 从星图中选出最亮的starNumMax颗星；
% 2. 让这starNumMax个对象星分别作为顶点星，遍历角距数据库，找到误差在±εd范围内的片段；
% 3. 对于每一个星组，标记每个片段中星等误差在±εb范围内的所有导航星并记录；
% 4. 通过确认匹配组中角距离关系，找出最大一致匹配组；
function [distMat, roughAtt] = AngMatch(AngLib, starGroup, starNumMax)
    starAngDistSet = zeros(starNumMax,starNumMax);
    for i = 1:starNumMax
        poleStar = [starGroup(i, 2), starGroup(i, 3), starGroup(i, 4)];
        for j = 1:starNumMax
            if i ~= j
                satellite = [starGroup(j, 2), starGroup(j, 3), starGroup(j, 4)];
                AngDist = acosd(dot(poleStar, satellite)/(vecnorm(poleStar)*vecnorm(satellite)));
                starAngDistSet(i, j) = AngDist;
                MatchGroup = find(AngLib(:, 1) < AngDist + func.matchconfig.errordist || AngLib(:, 1) > AngDist - func.matchconfig.errordist);

            else
                starAngDistSet(i, j) = 0;
            end
        end
    end
    distMat = starAngDistSet;
    roughAtt = [0,0];
end



function dcm = DualVecAttiCalc(libVec0, libVec1, senVec0, senVec1)
    % libVec0,1:    1*3 matrix, star vector in celestial coordinate system.
    % senVec0,1:    1*3 matrix, star vector in star sensor coordinate system.
    % dcm:          3*3 matrix, direction cosine matrix from J2000 to sensor.

    crossVec = cross(senVec0, senVec1);
    crossVec = crossVec / norm(crossVec);
    Ccm = [senVec0; crossVec; cross(senVec0, crossVec)];

    crossVec = cross(libVec0, libVec1);
    crossVec = crossVec / norm(crossVec);
    Ccr = [libVec0; crossVec; cross(libVec0, crossVec)];

    dcm = Ccm\Ccr;
end





