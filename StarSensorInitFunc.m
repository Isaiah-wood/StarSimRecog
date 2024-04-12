function func = StarSensorInitFunc
    func.WorkspacePath = strrep(fileparts(mfilename('fullpath')), '\', '/');
    func.HIPdatPath = [func.WorkspacePath, '/StarLib/hip_table.dat'];
    func.HIPcsvPath = [func.WorkspacePath, '/StarLib/hip_table.csv'];

    % func.StarLib = GetStarCoor(func.HIPdatPath,'xyz');
    % func.StarLib = GetStarMatrix(func.HIPdatPath,func.HIPcsvPath);
    func.StarLib = readmatrix(func.HIPcsvPath);

    func.GetStarCoor = @GetStarCoor;
    func.GetStarMatrix = @GetStarMatrix;

    func.Vec2Coo = @Vec2Coo;
    func.Coo2Vec = @Coo2Vec;
    func.Att2Dcm = @Att2Dcm;
    func.Dcm2Att = @Dcm2Att;
end

function pixCoo = Vec2Coo(sensorConf, dirVec)
    % Convert star point direction vector (in sensor coordinate system) to coordinate by pixel (in display coordinate system).
    % sensorConf:   sensor config.
    % dirVec:       3*n matrix, direction vectors in sensor coordinate system.
    % pixCoo:       2*n matrix, point coordinate in display coordinate system. Count from 1, the origin is upper left corner, by sequence of [row; col].

    % intrinsicParamMat = [sensorConf.f, 0, sensorConf.mainprow;...
    %                      0, sensorConf.f, sensorConf.mainpcol;...
    %                      0, 0, 1];
    intrinsicParamMat = [sensorConf.f, 0, sensorConf.width;...
                         0, sensorConf.f, sensorConf.height;...
                         0, 0, 1];
    % intrinsicParamMat = [sensorConf.f, 0, 0;...
    %                      0, sensorConf.f, 0;...
    %                      0, 0, 1];
    homoCoord = intrinsicParamMat * dirVec; % homogeneous coordinate
    pixCoo = homoCoord(1:2, :) ./ homoCoord(3, :) / sensorConf.pixelsize; % Point coordinate in image coordinate system. The origin is the main point. X axis is horizontal, Y axis is vertical.
end

function dirVec = Coo2Vec(sensorConf, pixCoo)
    % Convert star point coordinate by pixel (in display coordinate system) to direction vector (in sensor coordinate system).
    % sensorConf:   sensor config.
    % pixCoo:       2*n matrix, point coordinate in display coordinate system. Count from 1, the origin is upper left corner, by sequence of [row; col].
    % dirVec:       3*n matrix, direction vectors in sensor coordinate system.

    homoCoord = [pixCoo * sensorConf.pixelsize; ones(1, size(pixCooMillimeter_ImgCS, 2))]; % homogeneous coordinate
    intrinsicParamMatInv = [1 / sensorConf.f, 0, -sensorConf.mainprow / sensorConf.f;...
                            0, 1 / sensorConf.f, -sensorConf.mainpcol / sensorConf.f;...
                            0, 0, 1];
    dirVec = intrinsicParamMatInv * homoCoord;
    dirVec = dirVec ./ vecnorm(dirVec);
end

function dcm = Att2Dcm(att)
    % Convert star sensor attitude angle to direction cosine matrix.
    % Direction vectors 'vj' in J2000 could be converted to vectors 'vs' in sensor coordinate by 'vs = dcm * vj'.
    % att:  1*3 matrix, [ra dec roa], radius.
    % dcm:  3*3 matrix, direction cosine matrix from J2000 to sensor.

    cosRa = cos(att(1));
    sinRa = sin(att(1));
    cosDec = cos(att(2));

    sinDec = sin(att(2));
    cosRoa = cos(att(3));
    sinRoa = sin(att(3));
    dcm = [cosRoa, sinRoa, 0; -sinRoa, cosRoa, 0; 0, 0, 1] *...
          [1, 0, 0; 0, cosDec, sinDec; 0, -sinDec, cosDec] *...
          [cosRa, 0, sinRa; 0, 1, 0; -sinRa, 0, cosRa] *...
          [1, 0, 0; 0, 0, -1; 0, 1, 0] * [0, -1, 0; 1, 0, 0; 0, 0, 1];
end

% function dcm = Att2Dcm(att)
%     % Convert star sensor attitude angle to direction cosine matrix.
%     % Direction vectors 'vj' in J2000 could be converted to vectors 'vs' in sensor coordinate by 'vs = dcm * vj'.
%     % att:  1*3 matrix, [ra dec roa], radius.
%     % dcm:  3*3 matrix, direction cosine matrix from J2000 to sensor.
% 
%     % 原z轴向量  
%     Z_axis = [0; 0; 1];  
%     % 计算旋转轴（叉积）  
%     rotation_axis = cross(Z_axis, att);  
%     rotation_axis = rotation_axis / norm(rotation_axis);
% 
%     % 计算旋转角度（点积和反正切）  
%     theta = acos(dot(Z_axis , att));  
%     theta_squared = theta^2;  
% 
%     K = [0, -rotation_axis(3), rotation_axis(2);  
%          rotation_axis(3), 0, -rotation_axis(1);  
%          -rotation_axis(2), rotation_axis(1), 0];  
% 
%     matrix_R = eye(3)+sin(theta)*K+(1-cos(theta))*(K^2 +theta_squared*eye(3));          %罗德里格斯旋转公式
% 
%     if det(matrix_R) ~= 0
%         dcm = inv(matrix_R);
%     else
%         dcm = eye(3);
%     end
% 
% 
% end

function att = Dcm2Att(dcm)
    % Convert star sensor direction cosine matrix to attitude angle.
    % Direction vectors 'vj' in J2000 could be converted to vectors 'vs' in sensor coordinate by 'vs = dcm * vj'.
    % att:  1*3 matrix, [ra dec roa], radius.
    % dcm:  3*3 matrix, direction cosine matrix from J2000 to sensor.

    [ra, negDec, roa] = dcm2angle([0, 0, 1; -1, 0, 0; 0, -1, 0] * dcm, 'zyx'); % [0, 0, 1; -1, 0, 0; 0, -1, 0] = angle2dcm(pi/2, -pi/2, 0,'yzx') ^ -1
    att = [ra, -negDec, roa];
end






function StarsLib = GetStarCoor(starsLibPath,coorClass)
    starsTable = readtable(starsLibPath, "NumHeaderLines", 9, "VariableNamingRule", "preserve");
    starsTable = starsTable(1:end-1, :);
    Index = starsTable{:, 2};
    Direction = starsTable{:, 4};
    Vmag = starsTable{:, 5};
    % STAR_LENGTH = length(Index); 
    
    % 使用 cellfun 和 strsplit 分割每个 cell 中的字符串  
    splitStrings = cellfun(@strsplit, Direction, 'UniformOutput', false);  
    splitNumbers = cellfun(@str2double, splitStrings, 'UniformOutput', false);  
    coordinate = cell2mat(splitNumbers);  
    
    % 将右升交点从小时转换为度  RAhms
    longitude = (coordinate(:,1) + coordinate(:,2)/60 + coordinate(:,3)/3600) * (360/24); % 1小时 = 15度  
    % 将赤纬从度、分、秒转换为度  DEdms
    latitude = coordinate(:,4) + coordinate(:,5)/60 + coordinate(:,6)/3600;  
    
    if coorClass == "ll"
        StarsLib = [Index, longitude, latitude, Vmag];
    elseif coorClass == "xyz"
        cartesianX = cosd(longitude) .* cosd(latitude); 
        cartesianY = sind(longitude) .* cosd(latitude); 
        cartesianZ = sind(latitude);  
        CoorXYZ = [cartesianX, cartesianY, cartesianZ];
        % 使用norm函数计算每一行（即每个坐标点）的模  
        magnitudes = sqrt(sum(CoorXYZ.^2, 2));
          
        % 创建一个和coordinates同样大小的矩阵，用于存储单位化后的坐标  
        % normCoorXYZ = zeros(size(CoorXYZ));  
          
        % 单位化坐标：将每个坐标点除以其模  
        % 注意：我们需要避免除以零的情况，因此设置一个小阈值  
        tolerance = 1e-10; % 设置一个小的容忍值，用于避免除以零  
        magnitudes(magnitudes < tolerance) = tolerance; % 将接近零的模替换为容忍值  
        % 进行单位化操作  
        normCoorXYZ = bsxfun(@rdivide, CoorXYZ, magnitudes);  

        StarsLib = [Index, normCoorXYZ, Vmag];
    end

end



function StarsLib = GetStarMatrix(HIPdatPath,HIPcsvPath)
    starsTable = readtable(HIPdatPath, "NumHeaderLines", 9, "VariableNamingRule", "preserve");
    starsTable = starsTable(1:end-1, :);
    Index = starsTable{:, 2};
    Direction = starsTable{:, 4};
    Vmag = starsTable{:, 5};
    % STAR_LENGTH = length(Index); 
    
    % 使用 cellfun 和 strsplit 分割每个 cell 中的字符串  
    splitStrings = cellfun(@strsplit, Direction, 'UniformOutput', false);  
    splitNumbers = cellfun(@str2double, splitStrings, 'UniformOutput', false);  
    coordinate = cell2mat(splitNumbers);  
    
    % 将右升交点从小时转换为度  RAhms
    longitude = (coordinate(:,1) + coordinate(:,2)/60 + coordinate(:,3)/3600) * (360/24); % 1小时 = 15度  
    % 将赤纬从度、分、秒转换为度  DEdms
    latitude = coordinate(:,4) + coordinate(:,5)/60 + coordinate(:,6)/3600;  
    
    cartesianX = cosd(longitude) .* cosd(latitude); 
    cartesianY = sind(longitude) .* cosd(latitude); 
    cartesianZ = sind(latitude);  
    
    CoorXYZ = [cartesianX, cartesianY, cartesianZ];
    
% 之后考虑用CoorXYZ ./ vecnorm(CoorXYZ)实现
    magnitudes = sqrt(sum(CoorXYZ.^2, 2));
    tolerance = 1e-10; % 设置一个小的容忍值，用于避免除以零  
    magnitudes(magnitudes < tolerance) = tolerance; % 将接近零的模替换为容忍值   
    normCoorXYZ = bsxfun(@rdivide, CoorXYZ, magnitudes);  
      
    StarsLib = [Index, normCoorXYZ, Vmag, longitude, latitude];
    
    writematrix(StarsLib, HIPcsvPath);

end



