function h = steerVector3D(c, frequencies, mic_positions, theta_phi_pairs)
% ------ 计算导向矢量

%fprintf('\t------------------------------------------\n');
%fprintf('\tStart calculating 3D steering vector...\n');

% Setup scanning grid using grid_resolution and dimensions
N_mic = size(mic_positions, 1);
N_freqs = length(frequencies);

x_0 = mean(mic_positions, 1);

h = zeros(N_mic, size(theta_phi_pairs, 1), N_freqs);
% 注意，我们假设波长和频率都为1，即，c（波传播速度）也为1.
for K = 1:N_freqs
    k = 2*pi*frequencies(K)/c;
    for I = 1:N_mic
        p_vec = mic_positions(I,:); % - x_ac;  %z=0
        p_vec_norm = norm(p_vec);   %2-norm
        
        %r_vec = [cosd(phi_i)*sind(theta_i), sind(phi_i)*sind(theta_i),cosd(theta_i)];
        r_vec = [...
            cosd(theta_phi_pairs(:,2)).*sind(theta_phi_pairs(:,1)),...
            sind(theta_phi_pairs(:,2)).*sind(theta_phi_pairs(:,1)),...
            cosd(theta_phi_pairs(:,1))...
            ];
        if (p_vec_norm <= eps)
            cos_ang = 0;
        else
            %cos_ang = dot(r_vec,p_vec)/p_vec_norm;
            rp_dot = sum(r_vec .* p_vec, 2);
            cos_ang = rp_dot / p_vec_norm;
        end
        
        %a = exp(-1j * 2 * pi * f/c * p_vec_norm * cos_ang);
        h(I, :, K) = exp(-1i*k*p_vec_norm*cos_ang);
    end
end

%fprintf('\tFinished calculating steering vector!\n');
%fprintf('\t------------------------------------------\n');
end
