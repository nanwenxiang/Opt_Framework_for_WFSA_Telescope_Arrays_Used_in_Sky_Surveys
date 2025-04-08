function [snr,distance_fwhm,signalphoton,noisephoton] = ...
    calsnr11(sigmat,noisemat,psf_afterconv,pixelscale,camerapixcelsize)
    %It is a classical snr calcualtion function
    %Use the summation of all photons to divide the stadard
    %variation of noise
    distance_fwhm = calculateFWHM(psf_afterconv);
    neff_after_conv = 4*pi*((distance_fwhm*(pixelscale*1e-6)/camerapixcelsize)^2);
    snr = sum(sum(sigmat))/((sum(sum(sigmat))+neff_after_conv*mean(mean(noisemat)))^(0.5));
    signalphoton = sum(sum(sigmat));
    noisephoton = mean(mean(noisemat));
    % -------------------------------------------------------------------------
    function [max_matrix,distance] = add_half_high_circle(matrix_data,radio)
        [max_matrix,Index_x] = max(matrix_data);
        [max_matrix,Index_y] = max(max_matrix);
        [min_matrix,~] = min(matrix_data);
        [min_matrix,~] = min(min_matrix);
        Index_x = Index_x(Index_y);

        psf_copy = repmat(matrix_data,1,1);
        psf_copy(psf_copy<(max_matrix+min_matrix)*radio) = 0;
        psf_copy(psf_copy>=(max_matrix+min_matrix)*radio) = 1;

        psf_copy = edge(psf_copy);
        [m,n] = find(psf_copy == 1);
        distance = 0;
        for i = 1:length(m)
            distance_new = 2*((m(i)-Index_x)^2+(n(i)-Index_y)^2)^(1/2);
            if distance_new>distance
                distance = distance_new;
            end
        end
    end
% -------------------------------------------------------------------------
    function distance = scaleOfPsf(matrix_data,radio)
        [max_matrix,Index_x] = max(matrix_data);
        [max_matrix,Index_y] = max(max_matrix);
        Index_x = Index_x(Index_y);
        psf_copy = repmat(matrix_data,1,1);
        psf_copy(psf_copy<(max_matrix*radio)) = 0;
        psf_copy(psf_copy>=(max_matrix*radio)) = 1;
        psf_copy = edge(psf_copy);
        [m,n] = find(psf_copy == 1);
        distance = 0;
        for i = 1:length(m)
            distance_new = ((m(i)-Index_x)^2+(n(i)-Index_y)^2)^(1/2);
            if distance_new>distance
                distance = distance_new;
            end
        end
    end
    
    function distance=calculateFWHM(matrix_data)
        [FWHM_x, FWHM_y, fitback] = guass_fit(matrix_data);
        distance = mean([FWHM_x, FWHM_y], 'all');
    end

end
