function func = FuncStarImgSim
    func.InitConf = @InitConf;
    func.TakePhoto = @TakePhoto;
    func.StarLibInVision = @StarLibInVision;
    func.PrintPhoto = @PrintPhoto;
    func.SaveImgWithDir = @SaveImgWithDir;

    func.Vec2Coo = @Vec2Coo;
    func.Coo2Vec = @Coo2Vec;
    func.Att2Dcm = @Att2Dcm;
    func.Dcm2Att = @Dcm2Att;

    func.Sph2Vec = @Sph2Vec;
    func.Vec2Sph = @Vec2Sph;
end

function [cameraConf, noiseConf] = InitConf
    % 定义传感器参数的结构体
    % 参数说明：
    % cameraConf: 结构体，表示传感器参数.
    %       --ra        传感器的赤经(弧度).
    %       --dec       传感器的赤纬(弧度).
    %       --roa       传感器的滚转角(弧度).
    %       --f         光学系统的焦距(正数, mm).
    %       --fovradius 视场角半径(取对角线, 弧度).
    %       --pixelsize 像元尺寸(um).
    %       --blursigma 失焦模糊的高斯标准差.
    %       --constantc 常量C, 取决于系统的设计指标和曝光参数, 影响星点的绝对亮度.
    %       --height    图像像素高度.
    %       --width     图像像素宽度.
    %       --mainprow  像主点坐标row(像素).
    %       --mainpcol  像主点坐标col(像素).
    %       --channels  图像通道数(取1或3).
    cameraConf = struct( ...
        'ra', 90, ...
        'dec', 90, ...
        'roa', 0, ...
        'height', 4096, ...
        'width', 4096, ...
        'pixelsize', 2.74, ...
        'fovradius', 6, ...
        'magmax', 7.5, ...
        'maglimit', 8, ...
        'blursigma', 4, ...
        'constantc', 70000, ...
        'channels', 1 ...
    );
    cameraConf.ra       = rand() * 360;
    cameraConf.dec      = rand() * 180;
    cameraConf.roa      = rand() * 360;
    cameraConf.mainpcol = 0.5 * (cameraConf.height + 1);
    cameraConf.mainprow = 0.5 * (cameraConf.width + 1);
    cameraConf.f        = (cameraConf.height^2 + cameraConf.width^2)^0.5/2 * cameraConf.pixelsize / tand(cameraConf.fovradius);

    noiseConf = struct( ...
        'gauss', struct('enable', 0, 'mu', 0, 'sigma', 3), ...
        'fake', struct('enable', 0, 'marked', 0), ...
        'wake', struct('enable', 0, 'marked', 0) ...
    );
end

function VisibleStarList = StarLibInVision(StarLib, cameraConf)
    % StarLib  : Star library.
    % sensorRa : Right ascension of sensor.
    % sensorDec: Declination of sensor.
    % fovRadius: FoV radius of sensor.

    phi_Ocam = cameraConf.ra;
    theta_Ocam = cameraConf.dec;
    fovRadius = cameraConf.fovradius;
    magmax = cameraConf.magmax;

    libSize = size(StarLib, 1);
    VisibleStarList = zeros(100, 5);

    starCnt = 0;
    for i = 1:libSize
        phi_P = StarLib(i,6);
        theta_P = 90 - StarLib(i,7);
        if StarLib(i,5) < magmax && ...
            sind(theta_Ocam)*sind(theta_P)*cosd(phi_Ocam-phi_P) + ...
            cosd(theta_Ocam)*cosd(theta_P) > cosd(fovRadius)
            % 球面三角中的余弦定理：cos a = cos b * cos c + sin b sin c cos A
            starCnt = starCnt + 1;
            % 序号，坐标，星等
            VisibleStarList(starCnt, 1) = StarLib(i, 1);
            VisibleStarList(starCnt, 2:4) = StarLib(i, 2:4);
            VisibleStarList(starCnt, 5) = StarLib(i, 5);
        end
    end
    VisibleStarList(starCnt + 1:end, :) = [];
end


function [imgStarConf, VisibleStarListSorted] = TakePhoto(cameraConf, VisibleStarList)
    %   imgStarConf     结构体，每个字段均为，字段如下
    %       --type      星点是高斯分布:101
    %       --row,col   星点坐标(从1开始计)
    %       --max       星点高斯分布的最大值
    %       --sigma     星点高斯分布的标准差
    %                   当星点为高动态拖尾，参数...为...

    ra = cameraConf.ra;
    dec = cameraConf.dec;
    roa = cameraConf.roa;
    VisibleStarIdList = VisibleStarList(:, 1);
    VisibleStarVecList = VisibleStarList(:, 2:4);
    VisibleStarMagList = VisibleStarList(:, 5);


    dcm = Att2Dcm(deg2rad([ra, dec, roa]));
    sensorVisibleStarVecList = dcm * VisibleStarVecList'; % 方向矢量旋转至星敏坐标系
    starCoordPixel = Vec2Coo(cameraConf, sensorVisibleStarVecList)'; % 投影变换

    VisibleStarList(:, 6) = 1;
    VisibleStarList(:, 7:8) = round(starCoordPixel(:,1:2));

    % 剔除不在矩形图像范围内的星. Cull stars that are not within the bounds of the rectangular image.
    for starIdx = length(VisibleStarVecList):-1:1
        if starCoordPixel(starIdx, 1) < 1 || starCoordPixel(starIdx, 1) > cameraConf.width || starCoordPixel(starIdx, 2) < 1 || starCoordPixel(starIdx, 2) > cameraConf.height
            VisibleStarList(starIdx, 6) = 0;
            starCoordPixel(starIdx, :) = [];
            VisibleStarMagList(starIdx) = [];
            VisibleStarIdList(starIdx) = [];
        end
    end
    VisibleStarListSorted = sortrows(VisibleStarList, [6 5]);


    % Construct imgStarConf
    % type: 201 - gaussian target
    %       202 - high dynamic target
    imgStarConf.type = 101 + zeros(size(VisibleStarIdList));
    imgStarConf.sigma = cameraConf.blursigma + zeros(size(VisibleStarIdList));
    imgStarConf.id = VisibleStarIdList;
    imgStarConf.col = starCoordPixel(:, 1);
    imgStarConf.row = starCoordPixel(:, 2);
    imgStarConf.max = cameraConf.constantc ./ 2.512 .^ VisibleStarMagList;
end



function starImg = PrintPhoto(cameraConf, imgStarConf, noiseConf)
    %PrintPhoto 根据星点列表、噪声参数、背景参数生成星图
    %   加背景
    starImg = zeros(cameraConf.height, cameraConf.width);


    %   加星点
    starNum = size(imgStarConf.id, 1);
    starTypeList = imgStarConf.type;
    starMaxList = imgStarConf.max;
    starSigmaList = imgStarConf.sigma;
    starColList = imgStarConf.col;
    starRowList = imgStarConf.row;
    for starIdx = 1:starNum
        starMax = starMaxList(starIdx);
        starSigma = starSigmaList(starIdx);
        starCol = starColList(starIdx);
        starRow = starRowList(starIdx);
        % 确定星点窗口边界
        boxRadius = 1;
        while starMax * exp(- (boxRadius ^ 2 * 0.5 / starSigma ^ 2)) >= 0.5
            boxRadius = boxRadius + 1;
        end
        % 加星点
        if mod(starTypeList(starIdx), 100) == 1 % 高斯分布
            le = max(1, floor(starCol - boxRadius));
            ri = min(cameraConf.width, ceil(starCol + boxRadius));
            to = max(1, floor(starRow - boxRadius));
            bo = min(cameraConf.height, ceil(starRow + boxRadius));
            for row = to:bo
                for col = le:ri
                    dRow = row - starRow;
                    dCol = col - starCol;
                    starIntensity = starMax * exp(- (dCol ^ 2 + dRow ^ 2) * 0.5 / starSigma ^ 2);
                    starImg(row, col) = starImg(row, col) + starIntensity;
                end
            end
        end
    end

    if noiseConf.wake.enable
        for wakestarIdx = 1:7
            starMax = 600000 ./ 2.512 .^ (2*rand()+7);
            starSigma = 4;
            starCol = round(4096*rand());
            starRow = round(4096*rand());
            % 确定星点窗口边界
            boxRadius = 20;
            radius = 20;
            thickness = 2; % 方框粗细 
        
                % 确保星点坐标在图像范围内  
                starCol = max(1, min(starCol, cameraConf.width));  
                starRow = max(1, min(starRow, cameraConf.height));  
        
                % 计算方框边界  
                x1 = max(1, starCol - radius - thickness/2);  
                y1 = max(1, starRow - radius - thickness/2);  
                x2 = min(cameraConf.width, starCol + radius + thickness/2);  
                y2 = min(cameraConf.height, starRow + radius + thickness/2);  
        
                % 绘制方框  
                for row = y1:y2  
                    for col = x1:x2  
                        % 计算当前点到星点的距离  
                        d = sqrt((col - starCol).^2 + (row - starRow).^2);  
        
                        % 判断点是否在内框或外框上  
                        if (d >= radius - thickness/2) && (d <= radius + thickness/2)  
                            % 设置方框颜色为白色（或所需颜色）  
                            starImg(row, col) = starImg(row, col) + 255; % 假设图像是uint8类型  
                        end  
                    end  
                end  
        
        
            % 加星点
            le = max(1, floor(starCol - boxRadius));
            ri = min(cameraConf.width, ceil(starCol + boxRadius));
            to = max(1, floor(starRow - boxRadius));
            bo = min(cameraConf.height, ceil(starRow + boxRadius));
            for row = to:bo
                for col = le:ri
                    dRow = row - starRow;
                    dCol = col - starCol;
                    starIntensity = starMax * exp(- (dCol ^ 2 + dRow ^ 2) * 0.5 / starSigma ^ 2);
                    starImg(row, col) = starImg(row, col) + starIntensity;
                end
            end
        end
    end
    
    if noiseConf.fake.enable
        for fakestarIdx = 1:7
            starMax = 600000 ./ 2.512 .^ (2*rand()+5);
            starSigma = 4;
            starCol = round(4096*rand());
            starRow = round(4096*rand());
            % 确定星点窗口边界
            boxRadius = 20;
            radius = 48;
            thickness = 2; % 方框粗细 
        
                % 确保星点坐标在图像范围内  
                starCol = max(1, min(starCol, cameraConf.width));  
                starRow = max(1, min(starRow, cameraConf.height));  
        
                % 计算方框边界  
                x1 = max(1, starCol - radius - thickness/2);  
                y1 = max(1, starRow - radius - thickness/2);  
                x2 = min(cameraConf.width, starCol + radius + thickness/2);  
                y2 = min(cameraConf.height, starRow + radius + thickness/2);  
        
                % 绘制方框  
                for row = y1:y2  
                    for col = x1:x2  
                        % 计算当前点到星点的距离  
                        d = sqrt((col - starCol).^2 + (row - starRow).^2);  
        
                        % 判断点是否在内框或外框上  
                        if (d >= radius - thickness/2) && (d <= radius + thickness/2)  
                            % 设置方框颜色为白色（或所需颜色）  
                            starImg(row, col) = starImg(row, col) + 255; % 假设图像是uint8类型  
                        end  
                    end  
                end  
        
        
            % 加星点
            le = max(1, floor(starCol - boxRadius));
            ri = min(cameraConf.width, ceil(starCol + boxRadius));
            to = max(1, floor(starRow - boxRadius));
            bo = min(cameraConf.height, ceil(starRow + boxRadius));
            for row = to:bo
                for col = le:ri
                    dRow = row - starRow;
                    dCol = col - starCol;
                    starIntensity = starMax * exp(- (dCol ^ 2 + dRow ^ 2) * 0.5 / starSigma ^ 2);
                    starImg(row, col) = starImg(row, col) + starIntensity;
                end
            end
        end
    end

    %   加高斯噪声
    if noiseConf.gauss.enable
        starImg = starImg + normrnd(noiseConf.gauss.mu, noiseConf.gauss.sigma, cameraConf.height, cameraConf.width);
    end


    starImg = round(starImg);
    starImg = max(min(starImg, 255), 0);
    starImg = uint8(starImg);
end


function SaveImgWithDir(cameraConf, starImg, VisibleStarListSorted)
    ra  = cameraConf.ra;
    dec = cameraConf.dec;
    roa = cameraConf.roa;

    toolboxPath = strrep(fileparts(mfilename('fullpath')), '\', '/');
    fileName = ['ra',num2str(ra),'_dec',num2str(dec),'_roa',num2str(roa)];
    dirPath = [toolboxPath,'/Img/',fileName];
    imwrite(starImg,[dirPath,'.png'])
    writematrix(VisibleStarListSorted,[dirPath,'.csv'])
end



























function pixCoo = Vec2Coo(cameraConf, dirVec)
    % Convert star point direction vector (in sensor coordinate system) to coordinate by pixel (in display coordinate system).
    % cameraConf:   sensor config.
    % dirVec:       3*n matrix, direction vectors in sensor coordinate system.
    % pixCoo:       2*n matrix, point coordinate in display coordinate system. Count from 1, the origin is upper left corner, by sequence of [row; col].

    intrinsicParamMat = [cameraConf.f, 0, cameraConf.width/2*cameraConf.pixelsize;...
                        0, cameraConf.f, cameraConf.height/2*cameraConf.pixelsize;...
                        0, 0, 1];
    % intrinsicParamMat = [cameraConf.f, 0, 0;...
    %                      0, cameraConf.f, 0;...
    %                      0, 0, 1];
    homoCoord = intrinsicParamMat * dirVec; % homogeneous coordinate
    pixCoo = homoCoord(1:2, :) ./ homoCoord(3, :) / cameraConf.pixelsize; % Point coordinate in image coordinate system. The origin is the main point. X axis is horizontal, Y axis is vertical.
end

function dirVec = Coo2Vec(cameraConf, pixCoo)
    % Convert star point coordinate by pixel (in display coordinate system) to direction vector (in sensor coordinate system).
    % cameraConf:   sensor config.
    % pixCoo:       2*n matrix, point coordinate in display coordinate system. Count from 1, the origin is upper left corner, by sequence of [row; col].
    % dirVec:       3*n matrix, direction vectors in sensor coordinate system.

    homoCoord = [pixCoo * cameraConf.pixelsize; ones(1, size(pixCoo, 2))]; % homogeneous coordinate
    intrinsicParamMatInv = [1 / cameraConf.f, 0, -cameraConf.width/2*cameraConf.pixelsize / cameraConf.f;...
                            0, 1 / cameraConf.f, -cameraConf.height/2*cameraConf.pixelsize / cameraConf.f;...
                            0, 0, 1];
    dirVec = intrinsicParamMatInv * homoCoord;
    dirVec = dirVec ./ vecnorm(dirVec);
end

% function dirVec = Coo2Vec(cameraConf, pixCoo)
%     % Convert star point coordinate by pixel (in display coordinate system) to direction vector (in sensor coordinate system).
%     % cameraConf:   sensor config.
%     % pixCoo:       2*n matrix, point coordinate in display coordinate system. Count from 1, the origin is upper left corner, by sequence of [row; col].
%     % dirVec:       3*n matrix, direction vectors in sensor coordinate system.

%     pixCoo_ImgCS = pixCoo(2:-1:1, :) - [cameraConf.mainpcol; cameraConf.mainprow]; % Point coordinate in image coordinate system. The origin is the main point.
%     pixCooMillimeter_ImgCS = pixCoo_ImgCS * cameraConf.pixelsize / 1000; % The origin is the main point.
%     homoCoord = [pixCooMillimeter_ImgCS; ones(1, size(pixCooMillimeter_ImgCS, 2))]; % homogeneous coordinate
%     intrinsicParamMatInv = [1 / cameraConf.f, 0, 0; 0, 1 / cameraConf.f, 0; 0, 0, 1];
%     dirVec = intrinsicParamMatInv * homoCoord;
%     dirVec = dirVec ./ vecnorm(dirVec);
% end



function dcm = Att2Dcm(att)
    % Convert star sensor attitude angle to direction cosine matrix.
    % Direction vectors 'vj' in J2000 could be converted to vectors 'vs' in sensor coordinate by 'vs = dcm * vj'.
    % att:  1*3 matrix, [ra dec roa], radius.
    % dcm:  3*3 matrix, direction cosine matrix from J2000 to sensor.

    att = -att;
    % att(1) = -att(1);
    % att(2) = -att(2);
    % att(3) = -att(3);


    cosRa  = cos(att(1));
    sinRa  = sin(att(1));
    cosDec = cos(att(2));
    sinDec = sin(att(2));
    cosRoa = cos(att(3));
    sinRoa = sin(att(3));

    % dcm = [cosRa, sinRa, 0; -sinRa, cosRa, 0; 0, 0, 1] *...
    %       [1, 0, 0; 0, cosDec, sinDec; 0, -sinDec, cosDec] *...
    %       [cosRoa, sinRoa, 0; -sinRoa, cosRoa, 0; 0, 0, 1];

    dcm =   [cosRa, sinRa, 0; -sinRa, cosRa, 0; 0, 0, 1] *...
            [cosDec, 0, -sinDec; 0, 1, 0; sinDec, 0, cosDec] *...
            [cosRoa, sinRoa, 0; -sinRoa, cosRoa, 0; 0, 0, 1];

    dcm = dcm';

    % dcm = [1, 0, 0; 0, cosRa, sinRa; 0, -sinRa, cosRa] * [cosDec, 0, sinDec; 0, 1, 0; -sinDec, 0, cosDec];

    % dcm =   [1, 0, 0; 0, cosRa, sinRa; 0, -sinRa, cosRa] *...
    %     [cosDec, 0, sinDec; 0, 1, 0; -sinDec, 0, cosDec] *...
    %     [cosRoa, sinRoa, 0; -sinRoa, cosRoa, 0; 0, 0, 1];

end


function att = Dcm2Att(dcm)
    % Convert direction cosine matrix to star sensor attitude angle.
    % Direction vectors 'vj' in J2000 could be converted to vectors 'vs' in sensor coordinate by 'vs = dcm * vj'.
    % dcm:  3*3 matrix, direction cosine matrix from J2000 to sensor.
    % att:  1*3 matrix, [ra dec roa], radius.

    % Decompose DCM into rotations around each axis
    % Ra = atan2(dcm(2,3), dcm(1,3));
    % Dec = acos(dcm(3,3));
    % Roa = atan2(-dcm(3,2), dcm(3,1));

    Ra = atan2(dcm(3,2), dcm(3,1));
    Dec = acos(dcm(3,3));
    Roa = atan2(dcm(1,3), -dcm(2,3));

    if Ra < 0 
        Ra = 2*pi + Ra;
    end

    if Roa < 0 
        Roa = 2*pi + Roa;
    end

    % Return the negative angles as in the original Att2Dcm function
    att = [Ra, Dec, Roa];
end

function sph = Vec2Sph(vec)
    % Convert direction vector to spherical coordinate.
    % vec:   1*3 matrix, direction vector.
    % sph:   1*2 matrix, [ra dec], radius.

    sph(1) = atan2(vec(2), vec(1));
    sph(2) = atan2(vec(3), sqrt(vec(1)^2 + vec(2)^2));
end

function vec = Sph2Vec(sph)
    % Convert spherical coordinate to direction vector.
    % sph:   1*2 matrix, [ra dec], radius.
    % vec:   1*3 matrix, direction vector.

    vec(1) = cosd(sph(1)) * cosd(90 - sph(2));
    vec(2) = sind(sph(1)) * cosd(90 - sph(2));
    vec(3) = sind(90 - sph(2));

end