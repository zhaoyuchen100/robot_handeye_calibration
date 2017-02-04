
function hand_eye_calibration
robot_data = csvread('Robot_Points.csv',0,0);
robot_data = robot_data(:,1:6);
vicon_data = csvread('VICON_Points.csv',0,0);
count = 1;
for i = 1:size(robot_data,1)
    tmp_r = rotz(robot_data(i,6)*pi/180)*roty(robot_data(i,5)*pi/180)*rotx(robot_data(i,4)*pi/180);
    tmp_v = rotx(vicon_data(i,4))*roty(vicon_data(i,5))*rotz(vicon_data(i,6));
    A_(:,:,i) = [tmp_r,robot_data(i,1:3)';0 0 0 1];
    B_(:,:,i) = [tmp_v,vicon_data(i,1:3)';0 0 0 1];
end
% for i = 1:size(robot_data,1)-1
% %     for j = 1:size(robot_data,1)
%     B(:,:,count) = B_(:,:,i)\B_(:,:,i+1);
%     A(:,:,count) = A_(:,:,i)\A_(:,:,i+1);
%     count = count + 1;
% %     end
% end
for i = 1:size(robot_data,1)
    for j = 1:size(robot_data,1)
    B(:,:,count) = B_(:,:,i)\B_(:,:,j);
    A(:,:,count) = A_(:,:,i)\A_(:,:,j);
    count = count + 1;
    end
end

M = zeros(3,3);
C = [];
d = [];
ind = [];
valid_count = 0;
for i = 1:size(A,3)
    alpha(:,:,i) = calculate_log(A(1:3,1:3,i));
    beta(:,:,i) = calculate_log(B(1:3,1:3,i));
    alpha_vet(:,:,i) = extract_vect(alpha(:,:,i));
    beta_vet(:,:,i) = extract_vect(beta(:,:,i));
    if (abs(norm(alpha_vet(:,:,i))-norm(beta_vet(:,:,i)))<0.0001)
    t_ = beta(:,:,i)*alpha(:,:,i)';
    if (~isnan(t_(1)))
        ind = [ind,i];
        M = M + beta_vet(:,:,i)*alpha_vet(:,:,i)';
        valid_count = valid_count+1;
    end
    end
end
disp(['total count:=',num2str(count)]);
disp(['valid count:=', num2str(valid_count)]);
R_est = (M'*M)^(-1/2)*M'
A = A(:,:,ind);
B = B(:,:,ind);
for i =1:size(A,3)
    C = [C; eye(3)-A(1:3,1:3,i)];
    d = [d;A(1:3,4,i)-R_est*B(1:3,4,i)];
end
b_est = (C'*C)^(-1)*C'*d
norm(b_est)
X_est = [R_est,b_est;0 0 0 1]
end

function vect = extract_vect(mat)
vect = [mat(3,2),mat(1,3),mat(2,1)]';
end

function alpha = calculate_log(R)
trace_ = R(1,1)+R(2,2)+R(3,3);
psi = acos((trace_-1)/2);
alpha = psi*(R-R')/(2*sin(psi));
end


function [x] = importfile(workbookFile,sheetName,startRow,endRow)
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    startRow = 2;
    endRow = 10;
end

%% Import the data
data = xlsread(workbookFile, sheetName, sprintf('A%d:F%d',startRow(1),endRow(1)));
for block=2:length(startRow)
    tmpDataBlock = xlsread(workbookFile, sheetName, sprintf('A%d:F%d',startRow(block),endRow(block)));
    data = [data;tmpDataBlock]; %#ok<AGROW>
end

%% Allocate imported array to column variable names
TX = data(:,1);
TY = data(:,2);
TZ = data(:,3);
RX = data(:,4);
RY = data(:,5);
RZ = data(:,6);
x = [TX,TY,TZ,RX,RY,RZ];
end