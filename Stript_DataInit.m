WorkspacePath      = strrep(fileparts(mfilename('fullpath')), '\', '/');
HIPdatPath         = [WorkspacePath, '/StarLib/hip_table.dat'];
HIPcsvPath         = [WorkspacePath, '/StarLib/hip_table.csv'];
AngLibPath         = [WorkspacePath, '/StarLib/angle_database.csv'];
AngLibUnsortedPath = [WorkspacePath, '/StarLib/angle_database_unsorted.csv'];





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


