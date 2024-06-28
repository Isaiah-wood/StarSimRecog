WorkspacePath      = strrep(fileparts(mfilename('fullpath')), '\', '/');
HIPdatPath         = [WorkspacePath, '/StarLib/hip_table.dat'];
HIPcsvPath         = [WorkspacePath, '/StarLib/hip_table.csv'];
AngLibPath         = [WorkspacePath, '/StarLib/angle_database.csv'];
AngLibUnsortedPath = [WorkspacePath, '/StarLib/angle_database_unsorted.csv'];

StarLib = readmatrix(HIPcsvPath);


star1 = [0.0143057324404720	-0.0432309198475350	0.998962678776579];
star2 = [-0.0642880965762385	-0.00846363712656091	0.997895489259869];
Dist12 = acosd(dot(star1, star2)/(vecnorm(star1)*vecnorm(star2)));

star1fl = [0.305873261992882,0.708766519849817,-0.635681970750951];
star2fl = [0.298655347551951,0.704938718539656,-0.643324479933953];
Dist12fl = acosd(dot(star1fl, star2fl)/(vecnorm(star1fl)*vecnorm(star2fl)));

% disp([Dist12,Dist12fl,Dist12fl-Dist12])







% 创建球体  
[xSphere, ySphere, zSphere] = sphere(100); % 50 是网格的精细度  

% 绘制球体  
figure;  
surf(xSphere, ySphere, zSphere);  
axis equal; % 保持 x, y, z 轴的比例相等  
camlight; % 添加光源  
lighting gouraud; % 使用 Gouraud 阴影  

% 设置球体半径（这里使用 1 作为单位球）  
R = 1;  

% 定义三维向量并归一化  
v = star1fl; % 示例三维向量  
v_normalized = v / norm(v);  
v_scaled = R * v_normalized; % 缩放向量到球体半径  

% 在球上标记三维向量方向  
hold on;  
plot3(v_scaled(1), v_scaled(2), v_scaled(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);  
text(v_scaled(1), v_scaled(2), v_scaled(3), 'Vector Direction', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  

% 定义三维向量并归一化  
v = star2fl; % 示例三维向量  
v_normalized = v / norm(v);  
v_scaled = R * v_normalized; % 缩放向量到球体半径  

% 在球上标记三维向量方向  
hold on;  
plot3(v_scaled(1), v_scaled(2), v_scaled(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);  
text(v_scaled(1), v_scaled(2), v_scaled(3), 'Vector Direction', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  




% 定义三维向量并归一化  
v = star1; % 示例三维向量  
v_normalized = v / norm(v);  
v_scaled = R * v_normalized; % 缩放向量到球体半径  

% 在球上标记三维向量方向  
hold on;  
plot3(v_scaled(1), v_scaled(2), v_scaled(3), 'go', 'MarkerSize', 10, 'LineWidth', 2);  
text(v_scaled(1), v_scaled(2), v_scaled(3), 'Vector Direction', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  

% 定义三维向量并归一化  
v = star2; % 示例三维向量  
v_normalized = v / norm(v);  
v_scaled = R * v_normalized; % 缩放向量到球体半径  

% 在球上标记三维向量方向  
hold on;  
plot3(v_scaled(1), v_scaled(2), v_scaled(3), 'go', 'MarkerSize', 20, 'LineWidth', 2);  
text(v_scaled(1), v_scaled(2), v_scaled(3), 'Vector Direction', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  






























