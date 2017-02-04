
function hand_eye_calibration_demo
A1 = [-0.989992, -0.141120,0, 0;
        0.141120, -0.989992, 0, 0;
            0, 0, 1, 0;
                0, 0, 0, 1];
B1 = [-0.989992 -0.138307 0.028036 -26.9559;
         0.138307 -0.911449 0.387470 -96.1332;
             -0.028036 0.387470 0.921456 19.4872;
                 0 0 0 1];
A2 = [0.070737 0 0.997495 -400;
         0 1 0 0;
           -0.997495 0 0.070737 400;
             0 0 0 1];
 B2 = [0.070737 0.198172 0.997612 -309.543;
         -0.198172 0.963323 -0.180936 59.0244;
             -0.977612 -0.180936 0.107415 291.177;
                 0 0 0 1];
A(:,:,1) = [A1];
A(:,:,2) = [A2];
B(:,:,1) = [B1];
B(:,:,2) = [B2];

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
    if (abs(norm(alpha_vet(:,:,i))-norm(beta_vet(:,:,i)))<0.1)
%         i
%         abs(norm(alpha(:,:,i))-norm(beta(:,:,i)))
%   t_ = beta(:,:,i)*alpha(:,:,i)';
    if (~isnan(t_(1)))
        ind = [ind,i];
    %     beta(:,:,i)*alpha(:,:,i)'
        M = M + beta_vet(:,:,i)*alpha_vet(:,:,i)';
        valid_count = valid_count+1;
    end
    end
end
R_est = (M'*M)^(-1/2)*M'
% R_est = [1.0000         0         0
%          0    0.9801   -0.1983
%          0    0.1983    0.9801]
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

function alpha = calculate_log(R)
trace_ = R(1,1)+R(2,2)+R(3,3);
psi = acos((trace_-1)/2);
alpha = psi*(R-R')/(2*sin(psi));
end

function vect = extract_vect(mat)
vect = [mat(3,2),mat(1,3),mat(2,1)]';
end
