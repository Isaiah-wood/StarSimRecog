close all
clear

WorkspacePath      = strrep(fileparts(mfilename('fullpath')), '\', '/');
HIPdatPath         = [WorkspacePath, '/StarLib/hip_table.dat'];
HIPcsvPath         = [WorkspacePath, '/StarLib/hip_table.csv'];
AngLibPath         = [WorkspacePath, '/StarLib/angle_database.csv'];
AngLibUnsortedPath = [WorkspacePath, '/StarLib/angle_database_unsorted.csv'];

StarLib = readmatrix(HIPcsvPath);


cameraConf.maglimit = 7.5;
cameraConf.fovradius = 6;


AngLib = GetAngLib(AngLibPath, AngLibUnsortedPath, StarLib, cameraConf);



function AngLib = GetAngLib(AngLibPath, AngLibUnsortedPath, StarLib, cameraConf)
    StarIdList = StarLib(:, 1);
    StarVecList = StarLib(:, 2:4);
    StarMagList = StarLib(:, 5);
    StarLibMagMask = StarMagList <= cameraConf.maglimit;
    StarVecsubList = StarVecList(StarLibMagMask, :);
    StarIdsubList = StarIdList(StarLibMagMask, :);

    disp(size(StarVecsubList, 1))
    % writematrix(StarVecsubList, AngLibPath);
    LibIndex = 1;
    StarAngLib = zeros(10, 1);
    for i = 1:size(StarVecsubList, 1)
        Vec1 = StarVecsubList(i, :);
        for j = i+1:size(StarVecsubList, 1)
            Vec2 = StarVecsubList(j, :);
            AngDist = acosd(dot(Vec1, Vec2));
            if AngDist < cameraConf.fovradius * 2
            % if AngDist < cameraConf.fovradius / 3
                StarAngLib(LibIndex, 1) = AngDist;
                StarAngLib(LibIndex, 2) = StarIdsubList(i, :);
                StarAngLib(LibIndex, 3) = StarIdsubList(j, :);
                LibIndex = LibIndex + 1;
            end
        end
        disp(i)
        if mod(i, 1000) == 0 || i == size(StarVecsubList, 1)
            StarAngLib(LibIndex+1:end, :) = [];
            writematrix(StarAngLib, sprintf('StarLib/StarAngLib_%.5d.csv', i));
            LibIndex = 1;
            StarAngLib = zeros(100, 1);
        end
    end
    StarAngLib(LibIndex+1:end, :) = [];
    writematrix(StarAngLib, AngLibUnsortedPath);

    AngLib = sortrows(StarAngLib, 1, "ascend");
    writematrix(AngLib, AngLibPath);
end