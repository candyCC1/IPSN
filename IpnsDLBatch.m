function All_CBFs = IpnsDLBatch(A_steer, Rhat, ~, D)
    M = size(Rhat, 1);
    Dvect = diag(D);

    sigma_sq = Dvect(1);
    
    R_bar_hat = Rhat + sigma_sq * eye(M);

    %Tr_R_bar_hat = trace(R_bar_hat);
    %R_bar_hat = R_bar_hat / Tr_R_bar_hat;
    
    [V_bar_hat, D_bar_hat] = sorted_eig(R_bar_hat);
    SigmaBarHat = diag(D_bar_hat);
    
    k = trace(R_bar_hat) / M;
    
    J = size(A_steer, 2);       %the spectral point , or the grid count of direction of angle.
    
    if 0
        a1 = A_steer(:,1);

        Mmat = V_bar_hat;

        disp('norm(R_bar_hat - M*(D_bar_hat)* M^H) = ');
        disp(norm(R_bar_hat - Mmat*(D_bar_hat)*Mmat'));
        
        R_bar_hat_prime_1 = R_bar_hat + k * a1 * a1';
        %x = sqrt(k) * a1;
        %z = M'*x;
        x = a1;
        z = Mmat' * x;
        
        D = diag(z./abs(z));
        w = D' * z;
        %w = sqrt(k) * D' * z;
        
        normw = norm(w);
        w = w/normw;
        
        rho = k * normw * normw;

        disp('norm(R_bar_hat - M*D*(D_bar_hat)* D^H * M^H) = ');
        disp(norm(R_bar_hat - Mmat*D*(D_bar_hat)* D' * Mmat'));
        
        disp('norm(R_bar_hat_prime_1 - M*D*(D_bar_hat + rho * w*w^H)* D^H * M^H) = ');
        disp(norm(R_bar_hat_prime_1 - Mmat*D*(D_bar_hat + rho * w*w')* D' * Mmat'));

        disp('norm(k*x*x^H - M*D*(rho*w*w^H)* D^H * M^H) = ');
        disp(norm(k*x*x' - Mmat*D*(rho*w*w')* D' * Mmat'));
        
        %disp('norm(x*x^H - M*D*(w*w^H)* D^H * M^H) = ');
        %disp(norm(x*x' - M*D*(w*w')* D' * M'));

        %disp('norm(x*x^H - M^H*D*(w*w^H)* D^H * M) = ');
        %disp(norm(x*x' - M'*D*(w*w')* D' * M));
        
        [V_bar_hat_prime_1, D_bar_hat_prime_1] = sorted_eig(R_bar_hat_prime_1);
        SigmaBarHatPrime_1 = diag(D_bar_hat_prime_1);
        
        Sigma_N_diff_1 = SigmaBarHatPrime_1 - SigmaBarHat;
        disp('norm(Sigma_N_diff_1)');
        disp(norm(Sigma_N_diff_1));

        a2 = A_steer(:,10);
        
        R_bar_hat_prime_2 = R_bar_hat + k * a2 * a2';
        [V_bar_hat_prime_2, D_bar_hat_prime_2] = sorted_eig(R_bar_hat_prime_2);
        SigmaBarHatPrime_2 = diag(D_bar_hat_prime_2);
        
        Sigma_N_diff_2 = SigmaBarHatPrime_2 - SigmaBarHat;
        disp('norm(Sigma_N_diff_2)');
        disp(norm(Sigma_N_diff_2));
    end
    
    
    All_CBFs = zeros(J, M);
    for j = 1:J
        a = A_steer(:, j);
        R_bar_hat_prime = R_bar_hat + k * a * a';
        [V_bar_hat_prime, D_bar_hat_prime] = sorted_eig(R_bar_hat_prime);
        SigmaBarHatPrime = diag(D_bar_hat_prime);
        
        Sigma_N_diff = SigmaBarHatPrime - SigmaBarHat;

        for m = 1:M
            f_theta_prime = sum(Sigma_N_diff(end-m+1:end));
            All_CBFs(j, m) = 1/f_theta_prime;
        end
    end
    
    All_CBFs = real(All_CBFs);
end
