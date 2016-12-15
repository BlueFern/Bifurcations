% Plot J_KIR against Kp using v_i

% vi_max = vi_max(2:2:end);
% vi_min = vi_min(2:2:end);
% K_p = K_p(2:2:end);

% Extract different parts of graph
K_p_oscillations = smooth(K_p(2253:end), 0.1,'rloess').*1e3;
K_p_top = smooth(K_p(1:740), 0.1,'rloess').*1e3;
K_p_bottom = smooth(K_p(1110:2248), 0.1,'rloess').*1e3;

vi_max_smooth = smooth(vi_max(2253:end));
vi_min_smooth = smooth(vi_min(2253:end));
vi_top = smooth(vi_max(1:740));
vi_bottom = smooth(vi_max(1110:2248));

% Parameters
F_KIR = 7.5e2;
gamma_i = 1970;
z_1 = 4.5e-3;
z_2 = 112;
z_3 = 4.2e-4;
z_4 = 12.6;
z_5 = -7.4e-2;

% Channel resting potential
v_KIR_oscillations = z_1 .* K_p_oscillations - z_2;
v_KIR_top = z_1 .* K_p_top - z_2;
v_KIR_bottom = z_1 .* K_p_bottom - z_2;

% Channel Conductance
g_KIR_max = exp(z_5 * vi_max_smooth + z_3 * K_p_oscillations - z_4);
g_KIR_min = exp(z_5 * vi_min_smooth + z_3 * K_p_oscillations - z_4);
g_KIR_top = exp(z_5 * vi_top + z_3 * K_p_top - z_4);
g_KIR_bottom = exp(z_5 * vi_bottom + z_3 * K_p_bottom - z_4);

% KIR channel K+ flux
J_KIR_max = F_KIR .* g_KIR_max ./ gamma_i .* (vi_max_smooth - v_KIR_oscillations);
J_KIR_min = F_KIR .* g_KIR_min ./ gamma_i .* (vi_min_smooth - v_KIR_oscillations);
J_KIR_top = F_KIR .* g_KIR_top ./ gamma_i .* (vi_top - v_KIR_top);
J_KIR_bottom = F_KIR .* g_KIR_bottom ./ gamma_i .* (vi_bottom - v_KIR_bottom);

figure(1);
hold all
plot(K_p_oscillations./1e3, smooth(J_KIR_max, 0.1,'rloess'), 'r')
plot(K_p_oscillations./1e3, smooth(J_KIR_min, 0.1,'rloess'), 'r')
plot(K_p_top./1e3, smooth(J_KIR_top, 0.1,'rloess'), 'b')
plot(K_p_bottom./1e3, smooth(J_KIR_bottom, 0.1,'rloess'), 'b')
plot([0,30],[0,0],'-')
xlim([3 30])
xlabel('K_p (mM)')
title('J_{KIR} (\muM/s)')
