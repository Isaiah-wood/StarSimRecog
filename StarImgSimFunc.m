function func = StarImgSimFunc
    %STARIMGGENFUNC
    func.GetStar = @GetStar;
    func.StarLibInVision = @StarLibInVision;
    func.ExampleConf = @ExampleConf;
    func.ImgGen = @ImgGen;
end

function imgStarConf = GetStar(sif, sensorConf, starLib, StarLibInVision, constMag)
    % sensorConf:       结构体，表示传感器参数. Parameters of sensor.
    %       --ra        传感器的赤经(弧度). Right ascension of sensor (radian).
    %       --dec       传感器的赤纬(弧度). Declination of sensor (radian).
    %       --roa       传感器的滚转角(弧度). Roll angle of sensor (radian).
    %       --f         光学系统的焦距(正数, mm). Focal length of optical system (positive, mm).
    %       --fovradius 视场角半径(取对角线, 弧度). FoV radius of sensor (diagonal, radian).
    %       --pixelsize 像元尺寸(um). Pixel size (um).
    %       --blursigma 失焦模糊的高斯标准差. Gaussian standard deviation of out-of-focus blur.
    %       --constantc 常量C, 取决于系统的设计指标和曝光参数, 影响星点的绝对亮度. The constant C, depends on the design index and exposure parameters of the system, affects the absolute brightness of the star point.
    %       --height    图像像素高度. Height (pixels) of image.
    %       --width     图像像素宽度. Width (pixels) of image.
    %       --mainprow  像主点坐标row(像素, 从1计). Main point coordinate row (pixel, count from 1).
    %       --mainpcol  像主点坐标col(像素, 从1计). Main point coordinate col (pixel, count from 1).
    %       --channels  图像通道数(取1或3). Channels of image (1 or 3).
    %       注意: 焦距, 视场角半径, 像元尺寸和图像像素宽高之间存在制约关系, 知三求一.
    %       Note: there is a restrictive relationship between focal length, FoV radius, pixel size and image pixel width and height. Knowing the three can find the remaining one.
    % starLib:          Star library.
    % StarLibInVision:    Function handle to the function used to parse star library file.

    ra = sensorConf.ra;
    dec = sensorConf.dec;
    roa = sensorConf.roa;
    [VisibleStarVecList, VisibleStarMagList, VisibleStarIdList, VisibleStarSnList] = StarLibInVision(starLib, ra, dec, sensorConf.fovradius, sensorConf.magmax); % List of star direction vectors in celestial coordinate system.

    dcm = sif.Att2Dcm(deg2rad([ra, dec, roa]));
    sensorVisibleStarVecList = dcm * VisibleStarVecList; % 方向矢量旋转至星敏坐标系
    starCoordPixel = sif.Vec2Coo(sensorConf, sensorVisibleStarVecList); % 投影变换

    % 剔除不在矩形图像范围内的星. Cull stars that are not within the bounds of the rectangular image.
    % for starIdx = length(VisibleStarVecList):-1:1
    %     if starCoordPixel(1, starIdx) < 1 || starCoordPixel(1, starIdx) > sensorConf.height || starCoordPixel(2, starIdx) < 1 || starCoordPixel(2, starIdx) > sensorConf.width
    %         starCoordPixel(:, starIdx) = [];
    %         VisibleStarMagList(starIdx) = [];
    %         VisibleStarIdList(starIdx) = [];
    %         VisibleStarSnList(starIdx) = [];
    %     end
    % end
    starNum = length(VisibleStarVecList);


    % Cauculate star max brightness
    maxList = sensorConf.constantc ./ 2.512 .^ VisibleStarMagList;

    % Construct imgStarConf
    typeCell = 101 + zeros(starNum, 1);
    % type: 201 - gaussian target
    %       202 - high dynamic target
    imgStarConf.type = typeCell;
    imgStarConf.col = starCoordPixel(1, :)';
    imgStarConf.row = starCoordPixel(2, :)';
    imgStarConf.max = maxList';
    imgStarConf.sigma = sensorConf.blursigma + zeros(starNum, 1);
    imgStarConf.id = VisibleStarIdList';
    imgStarConf.sn = VisibleStarSnList';
end

function [VisibleStarVecList, VisibleStarMagList, VisibleStarIdList, VisibleStarSnList] = StarLibInVision(lib, sensorRa, sensorDec, fovRadius, magmax)
    % lib:          Star library.
    % sensorRa:     Right ascension of sensor.
    % sensorDec:    Declination of sensor.
    % fovRadius:    FoV radius of sensor.

    % VisibleStarVecList:      3xn matrix, list of direction vectors for possible stars in FoV, in celestial coordinate system.
    % VisibleStarMagList:    1xn matrix, list of magnitudes for possible stars in FoV.
    % VisibleStarIdList:       1xn matrix, list of identity for possible stars in FoV.
    % VisibleStarSnList:       1xn matrix, list of serial number for possible stars in FoV. Indicates the serial number of the star point in the star library file, count from zero.
    
    phi_Ocam = sensorRa;
    theta_Ocam = 90 - sensorDec;

    libSize = size(lib, 1);
    VisibleStarVecList = zeros(3, 100);
    VisibleStarMagList = zeros(1, 100);
    VisibleStarIdList = zeros(1, 100);
    VisibleStarSnList = zeros(1, 100);

    starCnt = 0;
    % 根据编号遍历数据  
    for i = 1:libSize
        phi_P = lib(i,6);
        theta_P = 90 - lib(i,7);
        if lib(i,5) < magmax && ...
        sind(theta_Ocam)*sind(theta_P)*cosd(phi_Ocam-phi_P) + ...
        cosd(theta_Ocam)*cosd(theta_P) > cosd(fovRadius)
        % 球面三角中的余弦定理：cos a = cos b * cos c + sin b sin c cos A
            starCnt = starCnt + 1;
            VisibleStarVecList(:, starCnt) = [lib(i, 2); lib(i, 3); lib(i, 4)];
            VisibleStarMagList(starCnt) = lib(i, 5);
            VisibleStarIdList(starCnt) = lib(i, 1);
            VisibleStarSnList(starCnt) = i - 1;
        end
    end

    VisibleStarVecList(:, starCnt + 1:end) = [];
    VisibleStarMagList(:, starCnt + 1:end) = [];
    VisibleStarIdList(:, starCnt + 1:end) = [];
    VisibleStarSnList(:, starCnt + 1:end) = [];
end

function [sensorConf, imgStarConf, imgBackgdConf, noiseConf] = ExampleConf
    sensorConf = struct( ...
        'ra', 0, ...
        'dec', 0, ...
        'roa', 0, ...
        'height', 4096, ...
        'width', 4096, ...
        'pixelsize', 2.74, ...
        'fovradius', 6, ...
        'magmax', 7.5, ...
        'maglimit', 8, ...
        'blursigma', 4, ...
        'constantc', 600000, ...
        'channels', 1 ...
    );
    sensorConf.ra = rand() * 2*pi;
    sensorConf.dec = rand() * pi;
    sensorConf.roa = rand() * 2*pi;
    sensorConf.f = (sensorConf.height^2 + sensorConf.width^2)^0.5/2 * sensorConf.pixelsize / tand(sensorConf.fovradius);
    sensorConf.mainpcol = 0.5 * (sensorConf.height + 1);
    sensorConf.mainprow = 0.5 * (sensorConf.width + 1);

    imgStarConf = struct( ...
        'type', [101; 101; 101; 101; 101; 101], ...
        'row', [354.4; 213.3; 10; 466.32; 135.6; 396.2], ...
        'col', [112.81; 432.1; 20; 608.66; 221.9; 503.4], ...
        'max', [120; 33; 80; 73; 60; 48], ...
        'sigma', [0.8; 0.8; 0.8; 0.8; 0.8; 0.8], ...
        'id', [0; 0; 0; 0; 0; 0], ...
        'sn', [0; 0; 0; 0; 0; 0] ...
    );
    imgBackgdConf = struct( ...
        'value', 30 ...
    );
    noiseConf = struct( ...
        'gauss', struct('enable', 1, 'mu', 0, 'sigma', 3) ...
        );
end

function [starImg, nbPixList] = ImgGen(sensorConf, imgStarConf, imgBackgdConf, noiseConf)
    %IMGGEN 根据星点列表、噪声参数、背景参数生成星图

    %   sensorConf:     结构体，表示传感器参数. Parameters of sensor.
    %       --ra        传感器的赤经(弧度). Right ascension of sensor (radian).
    %       --dec       传感器的赤纬(弧度). Declination of sensor (radian).
    %       --roa       传感器的滚转角(弧度). Roll angle of sensor (radian).
    %       --f         光学系统的焦距(正数, mm). Focal length of optical system (positive, mm).
    %       --fovradius 视场角半径(取对角线, 弧度). FoV radius of sensor (diagonal, radian).
    %       --pixelsize 像元尺寸(um). Pixel size (um).
    %       --blursigma 失焦模糊的高斯标准差. Gaussian standard deviation of out-of-focus blur.
    %       --constantc 常量C, 取决于系统的设计指标和曝光参数, 影响星点的绝对亮度. The constant C, depends on the design index and exposure parameters of the system, affects the absolute brightness of the star point.
    %       --height    图像像素高度. Height (pixels) of image.
    %       --width     图像像素宽度. Width (pixels) of image.
    %       --mainprow  像主点坐标row(像素). Main point coordinate row (pixel).
    %       --mainpcol  像主点坐标col(像素). Main point coordinate col (pixel).
    %       --channels  图像通道数(取1或3). Channels of image (1 or 3).
    %       注意: 焦距, 视场角半径, 像元尺寸和图像像素宽高之间存在制约关系, 知三求一.
    %       Note: there is a restrictive relationship between focal length, FoV radius, pixel size and image pixel width and height. Knowing the three can find the remaining one.

    %   imgStarConf     结构体，每个字段均为，字段如下
    %       --type      星点是高斯分布:101
    %       --row,col   星点坐标(从1开始计)
    %       --max       星点高斯分布的最大值
    %       --sigma     星点高斯分布的标准差
    %                   当星点为高动态拖尾，参数...为...

    %   backgdConf  结构体，表示背景参数，字段如下
    %       --value 背景灰度

    %   noiseConf           结构体，表示噪声参数，字段如下
    %       --gauss         高斯噪声参数结构体，字段如下
    %           --enable    是否含有高斯噪声
    %           --mu        高斯噪声均值
    %           --sigma     高斯噪声标准差

    %   返回值starImg是一个矩阵，指生成的星图，注意：为节省内存空间，无论sensorConf.channels取何值，starImg都为单通道图像
    % ----------------------------------------------------------------
    % 调用示例:
    % [starImg, nbPixList] = ImgGen(sensorConf, imgStarConf, imgBackgdConf, noiseConf);
    % imshow(starImg);
    % ----------------------------------------------------------------

    %   加背景
    starImg = zeros(sensorConf.height, sensorConf.width) + imgBackgdConf.value;

    %   初始化非背景点列表
    nbPixNum = 0;
    nbPixList = uint32(zeros(100, 3));

    %   加星点
    starNum = size(imgStarConf.row, 1);
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
            ri = min(sensorConf.width, ceil(starCol + boxRadius));
            to = max(1, floor(starRow - boxRadius));
            bo = min(sensorConf.height, ceil(starRow + boxRadius));
            for row = to:bo
                for col = le:ri
                    dRow = row - starRow;
                    dCol = col - starCol;
                    starIntensity = starMax * exp(- (dCol ^ 2 + dRow ^ 2) * 0.5 / starSigma ^ 2);
                    starImg(row, col) = starImg(row, col) + starIntensity;
                    % 记录非背景点列表
                    if starIntensity > noiseConf.gauss.sigma || starIntensity > 0.5 * starMax
                        nbPixNum = nbPixNum + 1;
                        if round(abs(dRow)) == 0 && round(abs(dCol)) == 0
                            nbPixList(nbPixNum, :) = [row, col, 1]; % 是星中心像素
                        else
                            nbPixList(nbPixNum, :) = [row, col, 2]; % 是星弥散斑像素
                        end
                    end
                end
            end
        end
    end

    %   加高斯噪声
    if noiseConf.gauss.enable
        starImg = starImg + normrnd(noiseConf.gauss.mu, noiseConf.gauss.sigma, sensorConf.height, sensorConf.width);
    end

    % 消去 nbPixList 中未用到的区域
    if size(nbPixList, 1) > nbPixNum
        nbPixList(nbPixNum + 1:size(nbPixList, 1), :) = [];
    end

    starImg = round(starImg);
    starImg = max(min(starImg, 255), 0);
    starImg = uint8(starImg);
end
