mSimulation;

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
v = VisibleStarList(:,2:4); % 示例三维向量  
v_normalized = v / norm(v);  
v_scaled = R * v_normalized; % 缩放向量到球体半径  
  
% 在球上标记三维向量方向  
hold on;  
plot3(v_scaled(1), v_scaled(2), v_scaled(3), 'ro', 'MarkerSize', 20, 'LineWidth', 2);  
text(v_scaled(1), v_scaled(2), v_scaled(3), 'Vector Direction', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  
  
ra = cameraConf.ra;
dec = cameraConf.dec;
roa = cameraConf.roa;
dcm = Fnsis.Att2Dcm(deg2rad([ra, dec, roa]));
sensorVisibleStarVecList = dcm * v';

% 定义三维向量并归一化  
v = sensorVisibleStarVecList(:,1); % 示例三维向量  
v_normalized = v / norm(v);  
v_scaled = R * v_normalized; % 缩放向量到球体半径  
  
% 在球上标记三维向量方向  
hold on;  
plot3(v_scaled(1), v_scaled(2), v_scaled(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);  
text(v_scaled(1), v_scaled(2), v_scaled(3), 'Vector Direction', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  


% % 定义经纬度坐标  
% lat = 90; % 纬度，单位为度  
% lon = 0; % 经度，单位为度  
% 
% % 将经纬度转换为三维笛卡尔坐标  
% lat_rad = deg2rad(lat);  
% lon_rad = deg2rad(lon);  
% x_latlon = R * cos(lat_rad) * cos(lon_rad);  
% y_latlon = R * cos(lat_rad) * sin(lon_rad);  
% z_latlon = R * sin(lat_rad);  
% 
% % 在球上标记经纬度对应的点  
% plot3(x_latlon, y_latlon, z_latlon, 'go', 'MarkerSize', 10, 'LineWidth', 2);  
% text(x_latlon, y_latlon, z_latlon, ['Lat: ', num2str(lat), ', Lon: ', num2str(lon)], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  


% 定义经纬度坐标  
lat = 90 - dec; % 纬度，单位为度  
lon = ra; % 经度，单位为度  

% 将经纬度转换为三维笛卡尔坐标  
lat_rad = deg2rad(lat);  
lon_rad = deg2rad(lon);  
x_latlon = R * cos(lat_rad) * cos(lon_rad);  
y_latlon = R * cos(lat_rad) * sin(lon_rad);  
z_latlon = R * sin(lat_rad);  

% 在球上标记经纬度对应的点  
plot3(x_latlon, y_latlon, z_latlon, 'go', 'MarkerSize', 20, 'LineWidth', 2);  
text(x_latlon, y_latlon, z_latlon, ['Lat: ', num2str(lat), ', Lon: ', num2str(lon)], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  

% 显示图形  
title('Sphere with Vector and Lat/Lon Markers');  
xlabel('X');  
ylabel('Y');  
zlabel('Z');  
grid on;  
hold off;