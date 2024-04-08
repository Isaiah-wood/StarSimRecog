function func = StarSensorBaseFunc
    func.toolboxPath = strrep(fileparts(mfilename('fullpath')), '\', '/');
    func.SAO60SSSLibPath = [func.toolboxPath, '/StarLib/sao60sss'];
    func.SAO60PAIRSLibPath = [func.toolboxPath, '/StarLib/PairsLib/sao60pairs_12'];
    func.starLib60 = load(func.SAO60SSSLibPath, '-ascii');
    func.pairLib60 = ReadPairsLib(func.SAO60PAIRSLibPath, 0.02 * pi / 180, 12 * pi / 180);

    func.Vec2Coo = @Vec2Coo;
    func.Coo2Vec = @Coo2Vec;
    func.Att2Dcm = @Att2Dcm;
    func.Dcm2Att = @Dcm2Att;
    func.CreatePairsLib = @CreatePairsLib;
    func.ReadPairsLib = @ReadPairsLib;
    func.PairsLibQuickIndex = @PairsLibQuickIndex;
end

function pixCoo = Vec2Coo(sensorConf, dirVec)
    % Convert star point direction vector (in sensor coordinate system) to coordinate by pixel (in display coordinate system).
    % sensorConf:   sensor config.
    % dirVec:       3*n matrix, direction vectors in sensor coordinate system.
    % pixCoo:       2*n matrix, point coordinate in display coordinate system. Count from 1, the origin is upper left corner, by sequence of [row; col].

    intrinsicParamMat = [sensorConf.f, 0, 0; 0, sensorConf.f, 0; 0, 0, 1];
    homoCoord = intrinsicParamMat * dirVec; % homogeneous coordinate
    pixCooMillimeter_ImgCS = homoCoord(1:2, :) ./ homoCoord(3, :); % The origin is the main point.
    pixCoo_ImgCS = pixCooMillimeter_ImgCS * 1000 / sensorConf.pixelsize; % Point coordinate in image coordinate system. The origin is the main point. X axis is horizontal, Y axis is vertical.
    pixCoo = pixCoo_ImgCS(2:-1:1, :) + [sensorConf.mainprow; sensorConf.mainpcol]; % Point coordinate in display coordinate system. Count from 1, the origin is upper left corner. X axis is vertical, Y axis is horizontal.
end

function dirVec = Coo2Vec(sensorConf, pixCoo)
    % Convert star point coordinate by pixel (in display coordinate system) to direction vector (in sensor coordinate system).
    % sensorConf:   sensor config.
    % pixCoo:       2*n matrix, point coordinate in display coordinate system. Count from 1, the origin is upper left corner, by sequence of [row; col].
    % dirVec:       3*n matrix, direction vectors in sensor coordinate system.

    pixCoo_ImgCS = pixCoo(2:-1:1, :) - [sensorConf.mainpcol; sensorConf.mainprow]; % Point coordinate in image coordinate system. The origin is the main point.
    pixCooMillimeter_ImgCS = pixCoo_ImgCS * sensorConf.pixelsize / 1000; % The origin is the main point.
    homoCoord = [pixCooMillimeter_ImgCS; ones(1, size(pixCooMillimeter_ImgCS, 2))]; % homogeneous coordinate
    intrinsicParamMatInv = [1 / sensorConf.f, 0, 0; 0, 1 / sensorConf.f, 0; 0, 0, 1];
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
    dcm = [cosRoa, sinRoa, 0; -sinRoa, cosRoa, 0; 0, 0, 1] * [1, 0, 0; 0, cosDec, sinDec; 0, -sinDec, cosDec] * [cosRa, 0, sinRa; 0, 1, 0; -sinRa, 0, cosRa] * [1, 0, 0; 0, 0, -1; 0, 1, 0] * [0, -1, 0; 1, 0, 0; 0, 0, 1];
end

function att = Dcm2Att(dcm)
    % Convert star sensor direction cosine matrix to attitude angle.
    % Direction vectors 'vj' in J2000 could be converted to vectors 'vs' in sensor coordinate by 'vs = dcm * vj'.
    % att:  1*3 matrix, [ra dec roa], radius.
    % dcm:  3*3 matrix, direction cosine matrix from J2000 to sensor.

    [ra, negDec, roa] = dcm2angle([0, 0, 1; -1, 0, 0; 0, -1, 0] * dcm, 'zyx'); % [0, 0, 1; -1, 0, 0; 0, -1, 0] = angle2dcm(pi/2, -pi/2, 0,'yzx') ^ -1
    att = [ra, -negDec, roa];
end

function CreatePairsLib(this, libName, fovDiagDist)
    % libName: new pairs library file name
    % fovDiagDist: Field of view diagonal distance (radian)
    lib = this.SAO60SSSLib;
    libSize = size(lib, 1);
    starVecList = [lib(:, 2), lib(:, 3), lib(:, 4)];
    magnitudeList = lib(:, 5);

    pairDist = zeros(200000, 1);
    pairId = zeros(200000, 2);
    pairNum = 0;

    pairsLibPath = sprintf("%s/StarLib/PairsLib/%s", this.sbf.toolboxPath, libName);

    fprintf(1, 'Scanning star %4d/%4d', 0, libSize);
    for i = 1:libSize
        fprintf(1, '\b\b\b\b\b\b\b\b\b%4d/%4d', i, libSize);
        vec_a = starVecList(i, :);
        for j = i + 1:libSize
            vec_b = starVecList(j, :);
            dist = acos(vec_a * vec_b' / sqrt((vec_a(1) * vec_a(1) + vec_a(2) * vec_a(2) + vec_a(3) * vec_a(3)) * (vec_b(1) * vec_b(1) + vec_b(2) * vec_b(2) + vec_b(3) * vec_b(3))));
            if dist <= fovDiagDist
                pairNum = pairNum + 1;
                pairDist(pairNum, 1) = dist;
                pairId(pairNum, 1) = i - 1;
                pairId(pairNum, 2) = j - 1;
            end
        end
    end
    fprintf(1, '\n');

    pairDist(pairNum + 1:end, :) = [];
    pairId(pairNum + 1:end, :) = [];

    [pairDist, sortedIdx] = sort(pairDist);
    pairId = pairId(sortedIdx, :);

    fid = fopen(pairsLibPath, 'w');
    for i = 1:pairNum
        fprintf(fid, '%1.8f %4d %4d\n', pairDist(i), pairId(i, 1), pairId(i, 2));
    end
    fclose(fid);
end

function pairsLib = ReadPairsLib(pairsLibPath, division, range)
    % division: represents the fineness of angular distance division (radian)
    % range：represents the maximum angular distance (radian)
    pairsList = load(pairsLibPath, '-ascii');
    grid = 0:division:range;
    idxList = zeros(size(grid, 2), 2);
    idxList(:, 2) = grid';
    pairIdx = 1;
    pairsNum = size(pairsList, 1);
    for gridIdx = 1:size(grid, 2)
        while pairIdx <= pairsNum
            if grid(gridIdx) <= pairsList(pairIdx, 1)
                break
            end
            pairIdx = pairIdx + 1;
        end
        idxList(gridIdx, 1) = pairIdx;
    end
    pairsLib = struct("pairs", pairsList, "index", idxList, "division", division, "range", range);
end

function [startIdx, endIdx] = PairsLibQuickIndex(pairsLib, dist, radius)
    % dist: angle distance (radian)
    % radius: represents the radius of the extended angular distance grid, default: 1

    if nargin < 3 || length(radius) == 0
        radius = 1;
    end

    idxList = pairsLib.index;
    idx = floor(dist / pairsLib.division) + 1;
    startIdx = idxList(min(size(idxList, 1), max(1, idx - radius)), 1);
    endIdx = idxList(min(size(idxList, 1), max(1, idx + 1 + radius)), 1) - 1;
end
