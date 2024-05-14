WorkspacePath      = strrep(fileparts(mfilename('fullpath')), '\', '/');
HIPdatPath         = [WorkspacePath, '/StarLib/hip_table.dat'];
HIPcsvPath         = [WorkspacePath, '/StarLib/hip_table.csv'];
AngLibPath         = [WorkspacePath, '/StarLib/angle_database.csv'];
AngLibUnsortedPath = [WorkspacePath, '/StarLib/angle_database_unsorted.csv'];

StarLib = readmatrix(HIPcsvPath);


VisibleStarList = readmatrix('VisibleStarList.csv')';

star1 = [0.034301808964359, -0.016085899558914, 0.999282057147607];
star2 = [0.054500623117466, 0.119014148736597, 0.991395639732348];
Dist12 = acosd(dot(star1, star2)/(vecnorm(star1)*vecnorm(star2)));

star1fl = [0.993293552468274, 0.110554179452966, 0.033848072772868];
star2fl = [0.999738556989801, -0.022785441479738, 0.001908225490767];
Dist12fl = acosd(dot(star1fl, star2fl)/(vecnorm(star1fl)*vecnorm(star2fl)));

disp([Dist12,Dist12fl,Dist12fl-Dist12])


VisibleStarListSorted = sortrows(VisibleStarList,5);


StarVecList = StarLib(:, 2:4);
StarMagList = StarLib(:, 5);
StarLibMagMask = StarMagList <= cameraConf.maglimit;
StarVecList = StarVecList(StarLibMagMask, : );




% img = starImg;
% imagesc(img); % 使用imagesc来显示图像，它会自动调整颜色映射以适应数据范围  
% colormap gray; % 设置为灰度颜色映射，因为你说这是一张黑白图片  
% axis image; % 设置为图像坐标轴，这样坐标轴就不会干扰图像显示  
% 
% text_x1 = 303;
% text_y1 = 1693;
% text(text_x1, text_y1, 'SN', 'FontSize', 12, 'Color', 'r');
% text_x2 = 4055.58642869495-10;
% text_y2 = 2262.85644968293;
% text(text_x2, text_y2, '2', 'Color', 'white', 'FontSize', 12);