function [g_out,q_out] = setup_matrix_abstract_v3()
syms eta_abs A_t1 A_t2 A_f1 A_f2 A_M1 A_M2 Q_2_inv eta_q P_r1 P_r2 th_c Q_M A_1 A_2 A_fl_2 A_fl_1

N = 6;
M = sym(zeros(N,N));

M(2,1) = eta_abs;
M(2,2) = eta_q*A_f1;
M(2,3) = P_r2*A_f2;
M(2,4) = P_r1*A_f2;

M(3,2) = 1/2*eta_q*(1-A_t1);
M(3,4) = P_r1*(1-A_t2);

M(4,2) = 1/2*eta_q*(1-A_t1);
M(4,3) = P_r2*(1-A_t2);

M(5,4) = 1-P_r1;% gain node

M(6,1) = 1-eta_abs; % loss node, gain+loss should equal one
M(6,2) = 1-eta_q+eta_q*A_M1;
M(6,3) = 1-P_r2+P_r2*A_M2;
M(6,4) = P_r1*A_M2;

M = subs(M, A_f1, A_t1 * 1 / (Q_2_inv + 1));
M = subs(M, A_f2, A_t2 * 1 / (Q_2_inv + 1));
M = subs(M, A_M1, A_t1 * Q_2_inv / (Q_2_inv + 1));
M = subs(M, A_M2, A_t2 * Q_2_inv / (Q_2_inv + 1));

M=simplify(M-eye(N));


M_inv = inv(M);
q=simplify(-M_inv(:,1));
G = simplify(-M_inv(N-1,1));

g_out = matlabFunction(G); % intensity at every node
q_out = matlabFunction(q);  % total efficiency of the system
% texm = latex(M) % for putting the solution into the report
% q_out(A,P_r1,P_r2,eta_q,eta_em,eta_abs) % use it like this
end