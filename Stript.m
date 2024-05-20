% WorkspacePath      = strrep(fileparts(mfilename('fullpath')), '\', '/');
% HIPdatPath         = [WorkspacePath, '/StarLib/hip_table.dat'];
% HIPcsvPath         = [WorkspacePath, '/StarLib/hip_table.csv'];
% AngLibPath         = [WorkspacePath, '/StarLib/angle_database.csv'];
% AngLibUnsortedPath = [WorkspacePath, '/StarLib/angle_database_unsorted.csv'];
% 
% StarLib = readmatrix(HIPcsvPath);
% 
% 
% VisibleStarList = readmatrix('VisibleStarList.csv')';
% 
% star1 = [0.034301808964359, -0.016085899558914, 0.999282057147607];
% star2 = [0.054500623117466, 0.119014148736597, 0.991395639732348];
% Dist12 = acosd(dot(star1, star2)/(vecnorm(star1)*vecnorm(star2)));
% 
% star1fl = [0.993293552468274, 0.110554179452966, 0.033848072772868];
% star2fl = [0.999738556989801, -0.022785441479738, 0.001908225490767];
% Dist12fl = acosd(dot(star1fl, star2fl)/(vecnorm(star1fl)*vecnorm(star2fl)));
% 
% disp([Dist12,Dist12fl,Dist12fl-Dist12])
% 
% 
% VisibleStarListSorted = sortrows(VisibleStarList,5);
% 
% 
% StarVecList = StarLib(:, 2:4);
% StarMagList = StarLib(:, 5);
% StarLibMagMask = StarMagList <= cameraConf.maglimit;
% StarVecList = StarVecList(StarLibMagMask, : );


