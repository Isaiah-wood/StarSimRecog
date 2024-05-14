function func = FuncStarMatch
    func.WorkspacePath = strrep(fileparts(mfilename('fullpath')), '\', '/');
    func.AngLibPath = [func.WorkspacePath, '/StarLib/angle_database.csv'];
    func.AngLib = readmatrix(func.AngLibPath);
    func.AngLibSize = size(func.AngLib, 1);

    func.matchconfig.starNumMax = 12;
    func.matchconfig.errordist = 0.03;
    func.matchconfig.errordbri = 0.05;

    func.AngLibSubErr = func.AngLib(:, 1) - func.matchconfig.errordist;
    func.AngLibAddErr = func.AngLib(:, 1) + func.matchconfig.errordist;


    func.AngMatch = @AngMatch;
    func.StarListFilter = @StarListFilter;
    func.DualVecAttiCalc = @DualVecAttiCalc;
end

function FilteredStarList = StarListFilter(starList, starNumMax)
    starBriList = starList.bri;
    sortedBriList = sort(starBriList, 'descend');
    briShel = sortedBriList(starNumMax+1);
    FilteredStarListMask = starBriList > briShel;
    FilteredStarList(:, 1) = starList.bri(FilteredStarListMask);
    FilteredStarList(:, 2) = starList.row(FilteredStarListMask);
    FilteredStarList(:, 3) = starList.col(FilteredStarListMask);
    FilteredStarList(:, 4) = starList.size(FilteredStarListMask);
    % FilteredStarList = sortrows(FilteredStarList, 1, "descend");


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




% function [distMat, roughAtt] = AngMatch(starGroup, starNumMax)
%     starAngDistSet = zeros(starNumMax,starNumMax);
%     for i = 1:starNumMax
%         poleStar = [starGroup(i, 2), starGroup(i, 3)];
%         for j = 1:starNumMax
%             if i ~= j
%                 satellite = [starGroup(j, 2), starGroup(j, 3)];
%                 AngDist = acosd(dot(poleStar, satellite)/(vecnorm(poleStar)*vecnorm(satellite)));
%                 starAngDistSet(i, j) = AngDist;
%             else
%                 starAngDistSet(i, j) = 0;
%             end
%         end
%     end
%     distMat = starAngDistSet;
%     roughAtt = [0,0];
% end



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





