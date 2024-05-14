function func = FuncStarRecog
    func.WorkspacePath = strrep(fileparts(mfilename('fullpath')), '\', '/');
    func.ImgDirPath = [func.WorkspacePath,'/Img/']; 

    func.ReadStarImg = @ReadStarImg;
    func.Binarization = @Binarization;
    func.CenterExtraction = @CenterExtraction;
    func.MarkPosition = @MarkPosition;
end

function [starImg, PriorAtt] = ReadStarImg(ImgDirPath, ImgName)
    % 读取指定目录下的星图图像，并提取图像名称中的先验姿态信息
    % ImgDirPath: 图像目录路径
    % ImgName: (可选) 图像文件名，如果提供则直接使用，否则随机选择
    
    ImgFiles = dir([ImgDirPath, '*.png']);
    
    % 检查是否存在PNG图像
    if isempty(ImgFiles)
        error('No PNG images found in the directory.');
    end

    % 如果未提供ImgName，随机选择一个图像
    if nargin < 2 || isempty(ImgName)
        fileIdx = randi(length(ImgFiles));
        fileName = ImgFiles(fileIdx).name;
    else
        fileName = ImgName; % 使用用户指定的图像文件名
    end

    StarImgPath = [ImgDirPath, fileName];

    % 提取图像名称中的先验姿态信息
    [~, name, ~] = fileparts(StarImgPath);
    pattern = '\d+(\.\d+)?';
    matches = regexp(name, pattern, 'match');
    PriorAtt = cellfun(@str2double, matches);

    % 读取图像
    starImg = imread(StarImgPath);

    % 显示结果（可选，根据需要保留或删除）
    disp(PriorAtt)
    % imshow(starImg)
end

function binaryImg = Binarization(srcImg, mode, arg1)
    % mode:
    % (default: MeanBased)
    % 'MeanBased' : arg1(default: 10) plus mean will be used as threshold of binarization.
    % 'FixedThres': arg1(default: 128) declares a fixed threshold of binarization.
    % 'TimesSigma': arg1(default: 3) declares that how many times sigma will be used as threshold.
    % 'RegionMeanBased' : binarization in a square region whose radius is arg2, arg1(default: 10) plus mean will be used as threshold of binarization.
    if nargin < 2
        mode = 'MeanBased';
    end
    if nargin < 3
        switch lower(mode)
            case 'fixedthres'
                arg1 = 128;
            case 'timessigma'
                arg1 = 3;
            case 'meanbased'
                arg1 = 10;
            case 'regionmeanbased'
                arg1 = 10;
            otherwise
                error(['There is no mode called ', mode, '.']);
                mode = 'MeanBased';
                arg1 = 10;
        end
    end

    switch lower(mode)
        case 'fixedthres'
            thres = arg1;
        case 'timessigma'
            imgMean = mean(mean(srcImg));
            imgStd = std2(srcImg);
            thres = imgMean + arg1 * imgStd;
        case 'meanbased'
            imgMean = mean(mean(srcImg));
            thres = imgMean + arg1;
        case 'regionmeanbased'
            arg1 = 10;
        otherwise
            error(['There is no mode called ', mode, '.']);
            imgMean = mean(mean(srcImg));
            thres = imgMean + arg1;
    end
    binaryImg = srcImg > thres;
end

function starList = CenterExtraction(srcImg, binaryImg, method, bgdThres, sizeThres)
    % method:
    % 'centroid'    质心法
    % 'sqweighted'  平方加权质心法
    % 'withthres'   带阈值的质心法，缺省时默认方法
    % 'fitting'     曲面拟合法，暂未实现

    % bgdThres:     当使用带阈值的质心法时，减去的背景值T，缺省时默认为输入图像的均值
    % sizeThres:    连通域的像素数量阈值, 尺寸小于此值的连通域将被剔除

    % starList:
    %   row:    1*n matrix
    %   col:    1*n matrix
    %   size:   1*n matrix
    %   bri:    1*n matrix

    imgMean = mean2(srcImg);
    if nargin < 3 || isempty(method)
        method = 'withthres';
    end
    if strcmp(method, 'withthres') && (nargin < 4 || isempty(bgdThres))
        bgdThres = imgMean;
    end
    if nargin < 5 || isempty(sizeThres)
        sizeThres = 1;
    end

    % 寻找所有连通区域
    connRegionStrc = bwconncomp(binaryImg, 4);
    connRegionList = connRegionStrc.PixelIdxList;
    % 尺寸阈值
    if sizeThres > 1
        starNum = 0;
        for starIdx = 1:length(connRegionList)
            if length(connRegionList{starIdx}) >= sizeThres
                starNum = starNum + 1;
                connRegionList{starNum} = connRegionList{starIdx};
            end
        end
        connRegionList(starNum + 1:end) = [];
    end
    starNum = length(connRegionList);
    % 宽高
    imgSize = size(srcImg);
    rows = imgSize(1);
    cols = imgSize(2);
    % 计算每个连通区域的重心
    % s = regionprops(cc, 'Centroid');
    % centroids = cat(1, s.Centroid);
    starSizeList = zeros(1, starNum);
    starCenRowList = zeros(1, starNum);
    starCenColList = zeros(1, starNum);
    starBrightnessList = zeros(1, starNum);
    if strcmp(method, 'centroid')
        for starIdx = 1:starNum
            % 对每一簇像素
            pixelList = connRegionList{starIdx};
            intensitySum = 0;
            cenRowSum = 0;
            cenColSum = 0;
            for pixelIdx = 1:length(pixelList)
                col = uint32(ceil(pixelList(pixelIdx) / rows));
                row = uint32(pixelList(pixelIdx) - (col - 1) * rows);
                intensitySum = intensitySum + double(srcImg(row, col));
                cenRowSum = cenRowSum + double(srcImg(row, col)) * double(row);
                cenColSum = cenColSum + double(srcImg(row, col)) * double(col);
            end
            cenRow = cenRowSum / intensitySum;
            cenCol = cenColSum / intensitySum;
            starCenRowList(starIdx) = cenRow;
            starCenColList(starIdx) = cenCol;
        end
    elseif strcmp(method, 'sqweighted')
        for starIdx = 1:starNum
            % 对每一簇像素
            pixelList = connRegionList{starIdx};
            intensitySum = 0;
            cenRowSum = 0;
            cenColSum = 0;
            for pixelIdx = 1:length(pixelList)
                col = uint32(ceil(pixelList(pixelIdx) / rows));
                row = uint32(pixelList(pixelIdx) - (col - 1) * rows);
                intensitySum = intensitySum + double(srcImg(row, col)) ^ 2;
                cenRowSum = cenRowSum + double(srcImg(row, col)) ^ 2 * double(row);
                cenColSum = cenColSum + double(srcImg(row, col)) ^ 2 * double(col);
            end
            cenRow = cenRowSum / intensitySum;
            cenCol = cenColSum / intensitySum;
            starCenRowList(starIdx) = cenRow;
            starCenColList(starIdx) = cenCol;
        end
    elseif strcmp(method, 'withthres')
        for starIdx = 1:starNum
            % 对每一簇像素
            pixelList = connRegionList{starIdx};
            intensitySum = 0;
            cenRowSum = 0;
            cenColSum = 0;
            for pixelIdx = 1:length(pixelList)
                col = uint32(ceil(pixelList(pixelIdx) / rows));
                row = uint32(pixelList(pixelIdx) - (col - 1) * rows);
                intensitySum = intensitySum + (double(srcImg(row, col)) - bgdThres);
                cenRowSum = cenRowSum + (double(srcImg(row, col)) - bgdThres) * double(row);
                cenColSum = cenColSum + (double(srcImg(row, col)) - bgdThres) * double(col);
            end
            cenRow = cenRowSum / intensitySum;
            cenCol = cenColSum / intensitySum;
            starCenRowList(starIdx) = cenRow;
            starCenColList(starIdx) = cenCol;
        end
    elseif strcmp(method, 'fitting')

    else
        ['Error: There are not a method called "', method, '".']
    end

    % 计算每颗星的大小及亮度（最亮四个像素的均值）
    for starIdx = 1:starNum
        % 对每一簇像素
        pixelList = connRegionList{starIdx};
        starSize = length(pixelList);
        starSizeList(starIdx) = starSize;
        if starSize >= 4
            sortedPix = sort(srcImg(pixelList), 'descend');
        else
            col = uint32(round(starCenColList(starIdx)));
            row = uint32(round(starCenRowList(starIdx)));
            win = srcImg(max(1, row - 1):min(rows, row + 1), max(1, col - 1):min(cols, col + 1));
            sortedPix = sort(reshape(win, [numel(win), 1]), 'descend');
        end
        starBrightnessList(starIdx) = mean(sortedPix(1:4)) - imgMean;
    end

    starList = struct('row', starCenRowList, 'col', starCenColList, 'size', starSizeList, 'bri', starBrightnessList);
end

function markImg = MarkPosition(starImg,starList)
    markImg = starImg;
    circle_mask = false(size(starImg));  
    % 遍历 row 和 col 数组，并在每个点处绘制一个半径为200的圆  
    for i = 1:length(starList.row)  
        % 计算当前圆心的坐标  
        cx = round(starList.col(i));  
        cy = round(starList.row(i));  
        [width,height] = size(starImg);
        radius = 40;
        thickness = 2; % 方框粗细 

        % 计算方框边界  
        x1 = max(1, cx - radius);  
        y1 = max(1, cy - radius);  
        x2 = min(width, cx + radius);  
        y2 = min(height, cy + radius);  
        
        % 遍历半径范围内的每个像素  
        for x = x1:x2  
            for y = y1:y2  
                % 计算当前像素到圆心的距离  
                d = sqrt((x - cx).^2 + (y - cy).^2);  
                % 如果像素在圆内或圆上（包括部分边界），则设置对应的 mask 为 true  
                if (d >= radius - thickness/2) && (d <= radius + thickness/2)  
                    circle_mask(y, x) = true;  
                end  
            end  
        end  
    end  


    markImg(circle_mask) = 255; % 设置为白色  

end