%%%%%%%%%%%%%%%%%%%%%%%%% 定义循环参数：拓扑flag和边长 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leaf = 1; sphere = 2; capsule = 3; ellipsoid = 4; concave = 5;
for flag = 1:5
    if flag == 1
        h0_values = [0.5,0.8,1,2,3];
    elseif flag == 2
		h0_values = [1,2,3,4,5];
	elseif flag == 3
		h0_values = [5,6,7,8.1,9];
    elseif flag == 4
        h0_values = [0.4,0.6,0.8,1,2];
	elseif flag == 5
		h0_values = [0.3,0.5,0.7,0.91,1.1]
    end
for h0 = h0_values

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 定义三维图形的几何参数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch flag
    case 1
        % 六瓣形（枫叶）
        A = [1/2, -sqrt(3)/2, 0; sqrt(3)/2, 1/2, 0; 0, 0, 1];%旋转矩阵绕 z 轴旋转60度
        B = [-1/2, -sqrt(3)/2, 0; sqrt(3)/2, -1/2, 0; 0, 0, 1];%旋转矩阵绕 z 轴旋转120度
        a = 48; b = 32; c = 8; %皮瓣椭球x轴、y轴、z轴半径长度
        d = 36; e =36; f = c; %胞体椭球半径长度
    case 2
        % 球形
        rs = 40;% 定义球半径
    case 3
        % 胶囊形
        r = 50; h = 1900; %半球半径和圆柱高度参数，胶囊纵向总长度为h+2*r，横向直径为2*r
    case 4
        % 椭球形
        ea = 80; eb = 48; ec = 8;  %皮瓣椭球x轴、y轴、z轴半径长度
    case 5
        % 双凹盘形
        D0 = 7.82; % 细胞直径
        a0 = 0.0518;a1 = 2.0026;a2 = -4.491;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 定义有向距离函数（SDF） %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch flag
    case 1
        % 六瓣形
        fd = @(p) dunion(dunion(dunion(sum((p).^2./repmat([a^2, b^2, c^2], size(p, 1), 1), 2) - 1, sum((p*A).^2./repmat([a^2, b^2, c^2], size(p, 1), 1), 2) - 1),...
									sum((p*B).^2./repmat([a^2, b^2, c^2], size(p, 1), 1), 2) - 1),sum((p).^2./repmat([d^2, d^2, c^2], size(p, 1), 1), 2) - 1);
		case 2
        % 球形
        fd = @(p) dsphere(p,0,0,0,rs);
    case 3
        % 胶囊形
        fd = @(p) p(:,1).^2 + p(:,2).^2 + (1/4) * (abs(p(:,3) - h/2) + abs(p(:,3) + h/2) - h).^2 - r^2;
    case 4
        % 椭球形
        fd = @(p) sum((p).^2./repmat([ea^2, eb^2, ec^2], size(p, 1), 1), 2)-1;    
    case 5
        % 双凹盘形
        fd = @(p)(D0 .* sqrt(1 - 4 .* (p(:,1).^2 + p(:,2).^2) / D0^2) .* (a0 + a1 .* (p(:,1).^2 + p(:,2).^2) / D0^2 + a2 .* (p(:,1).^2 + p(:,2).^2).^2 / D0^4)).^2 - p(:,3).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 定义网格生成参数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 定义边长函数
fh = @huniform;%均匀分布
% 定义盒子尺寸
switch flag
    case 1
				%六瓣形
        box0 = [-(max(a,b)+1),-(max(a,b)+1),-(c+1);(max(a,b)+1),(max(a,b)+1),(c+1)];
    case 2
				%球形
        box0 = [-rs,-rs,-rs; rs,rs,rs];
    case 3
				%胶囊形
        box0 = [-r,-r,-(2*r+h);r,r,(2*r+h)];
    case 4
				%椭球形
        box0 = [-2*ea,-2*eb,-2*ec;2*ea,2*eb,2*ec];
    case 5
				%双凹盘形
        box0 = [-D0,-D0,-D0;D0,D0,D0];
end
box = 1.1*box0;
% 定义初始边长
%h0 = 2;
% 生成网格并且计时
tic;
[p,t] = distmeshsurface(fd,fh,h0,box);
elapsedTime = toc;
% 寻找所有的边
e = [t(:, [1, 2]); t(:, [2, 3]); t(:, [3, 1])];
e = sort(e, 2);
e = unique(e, 'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 网格评估 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算网格单元数量
numElements = size(t, 1);
% 计算边的数量
numEdges = size(e, 1);
% 计算每个三角形的面积、法向量、角度
areas = zeros(size(t, 1), 1);
perimeters = zeros(size(t, 1), 1);
normals = zeros(size(t));
angles = zeros(size(t));
% 初始化总面积和总体积
areatotal = 0.0;
volumetotal = 0.0;
% 计算重心
gx = mean(p(:, 1));
gy = mean(p(:, 2));
gz = mean(p(:, 3));
for i = 1:size(t, 1)
    % 获取三角形的三个顶点的坐标
    v1 = p(t(i, 1), :);
    v2 = p(t(i, 2), :);
    v3 = p(t(i, 3), :);
    % 计算三角形两边的向量
    a = v2 - v1;
    b = v3 - v1;
    % a 和 b 的叉积给出了三角形的法向量
    normal = cross(a, b);
    normals(i, :) = normal;
    % 三角形的面积是法向量的模的一半
    areas(i) = norm(normal) / 2;
		% 计算每个三角形的周长
		perimeters(i) = norm(v2 - v1) + norm(v3 - v2) + norm(v1 - v3);
		% 处理角
		for j = 1:3
        % 获取三角形的三个顶点的坐标
        v1 = p(t(i, j), :);
        v2 = p(t(i, mod(j, 3) + 1), :);
        v3 = p(t(i, mod(j + 1, 3) + 1), :);
        % 计算 v1 处的角度
        angles(i, j) = acos(dot(v2 - v1, v3 - v1) / (norm(v2 - v1) * norm(v3 - v1)));
    end
		% 构造四维矩阵
    arr = [v1, gx;
           v2, gy;
           v3, gz;
           1, 1, 1, 1];
    % 计算并累加体积
    volumetotal = volumetotal + 1/6 * abs(det(arr));
end
% 计算每个三角形的内接圆半径和外接圆半径
inradii = areas ./ (0.5 * perimeters);
circumradii = zeros(size(t, 1), 1);
for i = 1:size(t, 1)
    p1 = p(t(i, 1), :);
    p2 = p(t(i, 2), :);
    p3 = p(t(i, 3), :);
    circumradii(i) = norm(p2 - p1) * norm(p3 - p2) * norm(p1 - p3) / (4 * areas(i));
end
% 计算每个三角形的质量
faceQualityRadius = inradii ./ circumradii;
faceQualityAreaPerimeter = 2 * sqrt(3) * areas ./ perimeters;
% 计算最小的面质量
minQualityRadius = min(faceQualityRadius);
minQualityAreaPerimeter = min(faceQualityAreaPerimeter);
% 计算面积的平均值
mean_area = mean(areas);
% 计算面积的标准差
std_area = std(areas);
% 变异系数（标准差除以平均值）是均匀性的一个度量
cv_area = std_area / mean_area;
% 转换为度
angles = angles * 180 / pi;
% 每个三角形的最小角度和最大角度是其质量的度量（最大角度越小，最小角度越大，质量越好）
min_angles = min(angles, [], 2);
max_angles = max(angles, [], 2);
% 计算最小角度和最大角度的平均值和标准差
mean_min_angle = mean(min_angles);
std_min_angle = std(min_angles);
mean_max_angle = mean(max_angles);
std_max_angle = std(max_angles);
% 计算每个顶点的连接数
vertexConnectivity = zeros(size(p, 1), 1);
for i = 1:size(t, 1)
    vertexConnectivity(t(i, :)) = vertexConnectivity(t(i, :)) + 1;
end
% 统计每个连接数的点的数量
edges = unique(vertexConnectivity);
counts = histc(vertexConnectivity, edges);
% 计算每个连接数的点的百分比
percentages = counts / sum(counts) * 100;
% 计算总表面积
areatotal = sum(areas);
% 计算所有边的长度
edgeLengths = sqrt(sum((p(e(:,2), :) - p(e(:,1), :)).^2, 2));
% 计算边长的最小值、最大值、平均值和标准差
minEdgeLength = min(edgeLengths);
maxEdgeLength = max(edgeLengths);
meanEdgeLength = mean(edgeLengths);
stdEdgeLength = std(edgeLengths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 输出配置 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 处理文件名
switch flag
    case 1
				modeltype = 'leaf';
    case 2
				modeltype = 'sphere';
    case 3
				modeltype = 'capsule';
    case 4
				modeltype = 'ellipsoid';
    case 5
				modeltype = 'concave';
end
filename = sprintf('result/%s_h0=%d', modeltype, h0);
if ~exist(filename, 'dir')
   mkdir(filename)
end
% 保存网格生成结果
dlmwrite(fullfile(filename, 'coord.txt'), p, 'delimiter', '\t');
dlmwrite(fullfile(filename, 'angle.txt'), t, 'delimiter', '\t');
dlmwrite(fullfile(filename, 'bond.txt'), e, 'delimiter', '\t');
saveas(gcf, fullfile(filename, 'mesh.png'));
% 保存盒子信息
fileID = fopen(fullfile(filename, 'box.txt'),'w');
fprintf(fileID, '%d %d\n', box(1,1),box(2,1));
fprintf(fileID, '%d %d\n', box(1,2),box(2,2));
fprintf(fileID, '%d %d\n', box(1,3),box(2,3));
fclose(fileID);
% 保存网格质量信息
fileID = fopen(fullfile(filename, 'statistic.txt'),'w');
fprintf(fileID, 'Number of Elements: %d\n', numElements);
fprintf(fileID, 'Number of edges: %d\n', numEdges);
fprintf(fileID, 'Minimum Mesh Quality by inradii / circumradii: %f\n', minQualityRadius);
fprintf(fileID, 'Minimum Mesh Quality by Area / Perimeter: %f\n', minQualityAreaPerimeter);
fprintf(fileID, 'Mesh Generation Time: %f\n', elapsedTime);
fprintf(fileID, 'Mean Min Angle: %f\n', mean_min_angle);
fprintf(fileID, 'Std Min Angle: %f\n', std_min_angle);
fprintf(fileID, 'Mean Max Angle: %f\n', mean_max_angle);
fprintf(fileID, 'Std Max Angle: %f\n', std_max_angle);
fprintf(fileID, 'Mean Area: %f\n', mean_area);
fprintf(fileID, 'Std Area: %f\n', std_area);
fprintf(fileID, 'CV Area: %f\n', cv_area);
fprintf(fileID, 'Surface area by triangle: %f\n', areatotal);
fprintf(fileID, 'Volume: %f\n', volumetotal);
fprintf(fileID, 'Minimum edge length: %f\n', minEdgeLength);
fprintf(fileID, 'Maximum edge length: %f\n', maxEdgeLength);
fprintf(fileID, 'Mean edge length: %f\n', meanEdgeLength);
fprintf(fileID, 'Standard deviation of edge length: %f\n', stdEdgeLength);
% 输出每个连接数的点的数量和百分比
for i = 1:length(counts)
    fprintf(fileID, 'Number of vertices with %d connections: %d\n', edges(i), counts(i));
    fprintf(fileID, 'Percentage of vertices with %d connections: %.2f%%\n', edges(i), percentages(i));
end
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end