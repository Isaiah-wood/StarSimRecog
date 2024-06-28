function func = FuncStarRecog
    func.WorkspacePath = strrep(fileparts(mfilename('fullpath')), '\', '/');
    func.ImgDirPath = [func.WorkspacePath,'/Img/']; 

    func.ReadStarImg = @ReadStarImg;
    func.Binarization = @Binarization;
    func.CenterExtraction = @CenterExtraction;
    func.MarkPosition = @MarkPosition;
end

function [starImg, priorAtt, starList] = ReadStarImg(ImgDirPath, FileName)
    % 读取指定目录下的星图图像，并提取图像名称中的先验姿态信息
    % ImgDirPath: 图像目录路径
    % FileName: (可选) 图像文件名，如果提供则直接使用，否则随机选择
    
    ImgFiles = dir([ImgDirPath, '*.png']);
    CsvFiles = dir([ImgDirPath, '*.csv']);
    
    % 检查是否存在PNG图像
    if isempty(ImgFiles)
        error('No PNG images found in the directory.');
    end
    if isempty(CsvFiles)
        error('No CSV images found in the directory.');
    end

    % 如果未提供FileName，随机选择一个图像
    if nargin < 2 || isempty(FileName)
        imgFileIdx = randi(length(ImgFiles));
        StarImgName = ImgFiles(imgFileIdx).name;
        csvFileIdx = randi(length(CsvFiles));
        StarCsvName = CsvFiles(csvFileIdx).name;
    else
        StarImgName = [FileName,'.png'];
        StarCsvName = [FileName,'.csv'];
    end

    StarImgPath = [ImgDirPath, StarImgName];
    StarCsvPath = [ImgDirPath, StarCsvName];

    % 提取图像名称中的先验姿态信息
    [~, name, ~] = fileparts(StarImgPath);
    pattern = '\d+(\.\d+)?';
    matches = regexp(name, pattern, 'match');
    priorAtt = cellfun(@str2double, matches);

    % 读取图像
    starImg = imread(StarImgPath);
    starList = readmatrix(StarCsvPath);
    starList = starList(starList(:,6)==1,:);

    % 显示结果（可选，根据需要保留或删除）
    disp(['读取星图的先验姿态:', num2str(priorAtt)])
    % imshow(starImg)
end

function binaryImg = Binarization(srcImg, mode, arg1)
    % 添加函数的文档字符串，说明函数的功能、输入输出参数等
    % binaryImg = Binarization(srcImg, mode, arg1) performs image binarization based on different modes.
    % @param srcImg Input image.
    % @param mode Binarization mode (default: 'MeanBased').
    % @param arg1 Argument depending on the mode.
    % @return binaryImg Binarized image.
    
    % 检查输入参数的有效性
    if ~isnumeric(srcImg) || isempty(srcImg)
        error('srcImg must be a numeric matrix.');
    end
    
    if nargin < 2
        mode = 'MeanBased';
    end
    
    % 设置默认的arg1值和检查mode的有效性合并进行
    defaultArg1 = struct('MeanBased', 10, 'FixedThres', 128, 'TimesSigma', 3);
    switch lower(mode)
        case {'meanbased', 'fixedthres', 'timessigma'}
            if nargin < 3
                arg1 = defaultArg1.(mode);
            end
        otherwise
            error(['There is no mode called ', mode, '.']);
    end
    
    % 根据不同的mode计算阈值
    switch lower(mode)
        case 'fixedthres'
            thres = arg1;
        case 'timessigma'
            times = arg1;
            imgMean = mean2(srcImg);
            imgStd = std2(srcImg);
            thres = imgMean + times * imgStd;
            disp(['thres:',num2str([imgMean, imgStd, thres])])
        case 'meanbased'
            offset = arg1;
            imgMean = mean2(srcImg);
            thres = imgMean + offset;
        otherwise
            error('Unknow mode.');
    end
    
    % 应用阈值生成二值图像
    binaryImg = srcImg > thres;
end

function ObservedStarList = CenterExtraction(srcImg, binaryImg, method, bgdThres, sizeThres)
    if nargin < 3 || isempty(method)
        method = 'withthres';
    end
    if strcmp(method, 'withthres') && (nargin < 4 || isempty(bgdThres))
        bgdThres = mean2(srcImg);
    end
    if nargin < 5 || isempty(sizeThres)
        sizeThres = 1;
    end

    % 计算图像均值（可能用于后续处理）
    imgMean = mean2(srcImg);
    
    % 使用bwconncomp寻找所有连通区域
    connRegionStrc = bwconncomp(binaryImg, 4);
    connRegionList = connRegionStrc.PixelIdxList;
    
    % 应用尺寸阈值
    if sizeThres > 1
        connRegionList = connRegionList(cellfun(@numel, connRegionList) >= sizeThres);
    end
    
    % 获取图像大小
    [rows,cols] = size(srcImg);
    srcImg = im2double(srcImg);
    % 初始化输出变量
    starNum = length(connRegionList);
    ObservedStarList = zeros(starNum, 4);

    for starIdx = 1:starNum
        pixelList = connRegionList{starIdx};
        % rowIndices = ceil(pixelList / rows);
        % colIndices = pixelList - (rowIndices - 1) * rows;
        [rowIndices, colIndices] = ind2sub([rows,cols],pixelList);
        
        % 根据方法选择不同的权重计算
        if strcmp(method, 'centroid')
            weight = srcImg(rowIndices, colIndices);
        elseif strcmp(method, 'sqweighted')
            weight = srcImg(rowIndices, colIndices).^2;
        elseif strcmp(method, 'withthres')
            weight = srcImg(rowIndices, colIndices) - bgdThres;
        end
        
        % 计算并更新结果
        cenRow = sum((rowIndices' * weight), "all") / sum(srcImg(rowIndices, colIndices), 'all');
        cenCol = sum(weight * colIndices, "all") / sum(srcImg(rowIndices, colIndices), 'all');
        starSize = numel(pixelList);
        starBrightness = mean2(srcImg(rowIndices, colIndices)) - imgMean;
        
        ObservedStarList(starIdx,:) = [cenRow, cenCol, starSize, starBrightness];
    end
end



% function ObservedStarList = CenterExtraction(srcImg, binaryImg, method, bgdThres, sizeThres)
%     % method:
%     % 'centroid'    质心法
%     % 'sqweighted'  平方加权质心法
%     % 'withthres'   带阈值的质心法，缺省时默认方法
%     % 'fitting'     曲面拟合法，暂未实现

%     % bgdThres:     当使用带阈值的质心法时，减去的背景值T，缺省时默认为输入图像的均值
%     % sizeThres:    连通域的像素数量阈值, 尺寸小于此值的连通域将被剔除

%     % ObservedStarList:
%     %   row:    1*n matrix
%     %   col:    1*n matrix
%     %   size:   1*n matrix
%     %   bri:    1*n matrix

%     imgMean = mean2(srcImg);
%     if nargin < 3 || isempty(method)
%         method = 'withthres';
%     end
%     if strcmp(method, 'withthres') && (nargin < 4 || isempty(bgdThres))
%         bgdThres = imgMean;
%     end
%     if nargin < 5 || isempty(sizeThres)
%         sizeThres = 1;
%     end

%     % 寻找所有连通区域
%     connRegionStrc = bwconncomp(binaryImg, 4);
%     connRegionList = connRegionStrc.PixelIdxList;
%     % 尺寸阈值
%     if sizeThres > 1
%         starNum = 0;
%         for starIdx = 1:length(connRegionList)
%             if length(connRegionList{starIdx}) >= sizeThres
%                 starNum = starNum + 1;
%                 connRegionList{starNum} = connRegionList{starIdx};
%             end
%         end
%         connRegionList(starNum + 1:end) = [];
%     end
%     starNum = length(connRegionList);
%     % 宽高
%     imgSize = size(srcImg);
%     rows = imgSize(1);
%     cols = imgSize(2);
%     % 计算每个连通区域的重心
%     % s = regionprops(cc, 'Centroid');
%     % centroids = cat(1, s.Centroid);
%     starSizeList       = zeros(starNum, 1);
%     starCenRowList     = zeros(starNum, 1);
%     starCenColList     = zeros(starNum, 1);
%     starBrightnessList = zeros(starNum, 1);
%     if strcmp(method, 'centroid')
%         for starIdx = 1:starNum
%             % 对每一簇像素
%             pixelList = connRegionList{starIdx};
%             intensitySum = 0;
%             cenRowSum = 0;
%             cenColSum = 0;
%             for pixelIdx = 1:length(pixelList)
%                 col = uint32(ceil(pixelList(pixelIdx) / rows));
%                 row = uint32(pixelList(pixelIdx) - (col - 1) * rows);
%                 intensitySum = intensitySum + double(srcImg(row, col));
%                 cenRowSum = cenRowSum + double(srcImg(row, col)) * double(row);
%                 cenColSum = cenColSum + double(srcImg(row, col)) * double(col);
%             end
%             cenRow = cenRowSum / intensitySum;
%             cenCol = cenColSum / intensitySum;
%             starCenRowList(starIdx) = cenRow;
%             starCenColList(starIdx) = cenCol;
%         end
%     elseif strcmp(method, 'sqweighted')
%         for starIdx = 1:starNum
%             % 对每一簇像素
%             pixelList = connRegionList{starIdx};
%             intensitySum = 0;
%             cenRowSum = 0;
%             cenColSum = 0;
%             for pixelIdx = 1:length(pixelList)
%                 col = uint32(ceil(pixelList(pixelIdx) / rows));
%                 row = uint32(pixelList(pixelIdx) - (col - 1) * rows);
%                 intensitySum = intensitySum + double(srcImg(row, col)) ^ 2;
%                 cenRowSum = cenRowSum + double(srcImg(row, col)) ^ 2 * double(row);
%                 cenColSum = cenColSum + double(srcImg(row, col)) ^ 2 * double(col);
%             end
%             cenRow = cenRowSum / intensitySum;
%             cenCol = cenColSum / intensitySum;
%             starCenRowList(starIdx) = cenRow;
%             starCenColList(starIdx) = cenCol;
%         end
%     elseif strcmp(method, 'withthres')
%         for starIdx = 1:starNum
%             % 对每一簇像素
%             pixelList = connRegionList{starIdx};
%             intensitySum = 0;
%             cenRowSum = 0;
%             cenColSum = 0;
%             for pixelIdx = 1:length(pixelList)
%                 col = uint32(ceil(pixelList(pixelIdx) / rows));
%                 row = uint32(pixelList(pixelIdx) - (col - 1) * rows);
%                 intensitySum = intensitySum + (double(srcImg(row, col)) - bgdThres);
%                 cenRowSum = cenRowSum + (double(srcImg(row, col)) - bgdThres) * double(row);
%                 cenColSum = cenColSum + (double(srcImg(row, col)) - bgdThres) * double(col);
%             end
%             cenRow = cenRowSum / intensitySum;
%             cenCol = cenColSum / intensitySum;
%             starCenRowList(starIdx) = cenRow;
%             starCenColList(starIdx) = cenCol;
%         end
%     elseif strcmp(method, 'fitting')

%     else
%         ['Error: There are not a method called "', method, '".']
%     end

%     % 计算每颗星的大小及亮度（最亮四个像素的均值）
%     for starIdx = 1:starNum
%         % 对每一簇像素
%         pixelList = connRegionList{starIdx};
%         starSize = length(pixelList);
%         starSizeList(starIdx) = starSize;
%         if starSize >= 4
%             sortedPix = sort(srcImg(pixelList), 'descend');
%         else
%             col = uint32(round(starCenColList(starIdx)));
%             row = uint32(round(starCenRowList(starIdx)));
%             win = srcImg(max(1, row - 1):min(rows, row + 1), max(1, col - 1):min(cols, col + 1));
%             sortedPix = sort(reshape(win, [numel(win), 1]), 'descend');
%         end
%         starBrightnessList(starIdx) = mean(sortedPix(1:4)) - imgMean;
%     end

%     % ObservedStarList = struct('row', starCenRowList, 'col', starCenColList, 'size', starSizeList, 'bri', starBrightnessList);
%     ObservedStarList = [starCenRowList, starCenColList, starSizeList, starBrightnessList];
% end




















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