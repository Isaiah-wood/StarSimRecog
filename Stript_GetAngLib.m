close all
clear
format longG

fovRadius = 6;

WorkspacePath    = strrep(fileparts(mfilename('fullpath')), '\', '/');
HIPdatPath       = [WorkspacePath, '/StarLib/hip_table.dat'];
HIPcsvPath       = [WorkspacePath, '/StarLib/hip_table.csv'];
AngLib_Mag8_Path       = [WorkspacePath, '/StarLib/AngLib_mag8.csv'];
AngLib_Mag8_SortedPath = [WorkspacePath, '/StarLib/AngLib_mag8_sorted.csv'];
AngLib_Mag7_5_Path       = [WorkspacePath, '/StarLib/AngLib_mag7_5.csv'];
AngLib_Mag7_5_SortedPath = [WorkspacePath, '/StarLib/AngLib_mag7_5sorted.csv'];


if isempty(dir(HIPcsvPath))
    StarsLib = dat2csv(HIPdatPath, HIPcsvPath);
else
    StarsLib = readmatrix(HIPcsvPath);
end


SubStarsLib_Mag7_5 = StarsLib(StarsLib(:, 5) <= 7.5, 1:4);
clear StarsLib;
AngLib_Mag7_5 = GetAngLib(AngLib_Mag7_5_Path, AngLib_Mag7_5_SortedPath, SubStarsLib_Mag7_5, fovRadius);



SubStarsLib_Mag8   = StarsLib(StarsLib(:, 5) <= 8, 1:4);
clear StarsLib;
AngLib_Mag8 = GetAngLib(AngLib_Mag8_Path, AngLib_Mag8_SortedPath, SubStarsLib_Mag8, fovRadius);








function AngLib = GetAngLib(AngLibPath, AngLibSortedPath, StarLib, fovRadius)
    disp(size(StarLib, 1))
    StarLibNum = size(StarLib, 1);
    StarIdList = StarLib(:, 1);
    StarVecList = StarLib(:, 2:4);
    
    LibIndex = 1;
    AngDistThreshold = cosd(fovRadius * 2);
    StarAngLib = zeros(10, 1);
    for i = 1:StarLibNum
        Vec1 = StarVecList(i, :);
        for j = i+1:StarLibNum
            Vec2 = StarVecList(j, :);
            AngDist = dot(Vec1, Vec2);
            if AngDist > AngDistThreshold
                StarAngLib(LibIndex, 1:3) = [AngDist, StarIdList(i, :), StarIdList(j, :)];
                LibIndex = LibIndex + 1;
            end
        end
        disp(i)
        if mod(i, 1000) == 0 || i == StarLibNum
            StarAngLib(LibIndex+1:end, :) = [];
            writematrix(StarAngLib, sprintf('StarLib/temp/StarAngLib_%.5d.csv', i));
            LibIndex = 1;
            StarAngLib = zeros(10, 1);
        end
    end
    StarAngLib(LibIndex+1:end, :) = [];

    ImgFiles = dir('StarLib/temp/StarAngLib_*.csv');
    if isempty(ImgFiles)
        disp('No PNG images found in the directory.')
    else
        StarAngLibAll = [];
        for i = 1:length(ImgFiles)
            filePath = fullfile('StarLib/temp/', ImgFiles(i).name);
            try
                StarAngLib = readmatrix(filePath);
            catch ME
                disp(['Error reading ' filePath ': ' ME.message]);
                continue;
            end
            
            % 删除全为0的行
            zeroRows = all(StarAngLib == 0, 2);  
            StarAngLib(zeroRows, :) = [];
            
            % 将处理后的数据存入cell数组
            % StarAngLibAll{end+1} = StarAngLib;
            StarAngLibAll = [StarAngLibAll; StarAngLib];
            delete(filePath); 
        end
        StarAngLib = StarAngLibAll;
    end

    writematrix(StarAngLib, AngLibSortedPath);
    AngLib = sortrows(StarAngLib, 1, "descend");
    writematrix(AngLib, AngLibPath);
end