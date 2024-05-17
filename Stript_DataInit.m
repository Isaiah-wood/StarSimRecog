%% 脚本————数据预处理
% 从HIP星库文件hip_table.dat中读取数据；
% 生成坐标子表，Vec表示三维向量坐标，Sph表示赤经赤纬
% 生成星等子表，从8/7.5星等处截取子表




WorkspacePath = strrep(fileparts(mfilename('fullpath')), '\', '/');
HIPdatPath    = [WorkspacePath, '/StarLib/hip_table.dat'];
HIPcsvPath    = [WorkspacePath, '/StarLib/hip_table.csv'];
VecLibPath    = [WorkspacePath, '/StarLib/starlib_vec.csv'];
SphLibPath    = [WorkspacePath, '/StarLib/starlib_sph.csv'];

SubLib_Mag8_Path   = [WorkspacePath, '/StarLib/sublib_mag8.csv'];
SubLib_Mag7_5_Path = [WorkspacePath, '/StarLib/sublib_mag7_5.csv'];


StarsLib    = dat2csv(HIPdatPath, HIPcsvPath);
StarsVecLib = StarsLib(:, 1:5);
StarsSphLib = [StarsLib(:, 1), StarsLib(:, 6:7), StarsLib(:, 5)];
writematrix(StarsVecLib, VecLibPath);
writematrix(StarsSphLib, SphLibPath);

SubStarsLib_Mag8   = StarsLib(StarsLib(:, 5) <= 8, :);
SubStarsLib_Mag7_5 = StarsLib(StarsLib(:, 5) <= 7.5, :);
writematrix(SubStarsLib_Mag8, SubLib_Mag8_Path);
writematrix(SubStarsLib_Mag7_5, SubLib_Mag7_5_Path);






%% 函数 dat2csv
function StarsLib = dat2csv(datPath, csvPath)

    starsTable = readtable(datPath,"NumHeaderLines",9,"VariableNamingRule","preserve");
    starsTable = starsTable(1:end-1, :);
    Index = starsTable{:, 2};
    Direction = starsTable{:, 4};
    Vmag = starsTable{:, 5};

    % 使用 cellfun 和 strsplit 分割每个 cell 中的字符串  
    splitStrings = cellfun(@strsplit, Direction, 'UniformOutput', false);
    splitNumbers = cellfun(@str2double, splitStrings, 'UniformOutput', false);
    coordinate   = cell2mat(splitNumbers);

    % 将右升交点从小时转换为度  RAhms
    longitude = (coordinate(:,1) + coordinate(:,2)/60 + coordinate(:,3)/3600) * (360/24); % 1小时 = 15度  
    % 将赤纬从度、分、秒转换为度  DEdms
    latitude = coordinate(:,4) + coordinate(:,5)/60 + coordinate(:,6)/3600;  

    cartesianX = cosd(longitude) .* cosd(latitude); 
    cartesianY = sind(longitude) .* cosd(latitude); 
    cartesianZ = sind(latitude);  

    CoorXYZ = [cartesianX, cartesianY, cartesianZ];

    % 使用norm函数计算每一行（即每个坐标点）的模  
    magnitudes = sqrt(sum(CoorXYZ.^2, 2));
    tolerance = 1e-10; % 设置一个小的容忍值，用于避免除以零  
    magnitudes(magnitudes < tolerance) = tolerance; % 将接近零的模替换为容忍值  

    % 创建一个和coordinates同样大小的矩阵，用于存储单位化后的坐标  
    % normCoorXYZ = zeros(size(CoorXYZ));  
    
    % 进行单位化操作  
    normCoorXYZ = bsxfun(@rdivide, CoorXYZ, magnitudes);  
    
    StarsLib = [Index, normCoorXYZ, Vmag, longitude, latitude];

    writematrix(StarsLib, csvPath);

end


