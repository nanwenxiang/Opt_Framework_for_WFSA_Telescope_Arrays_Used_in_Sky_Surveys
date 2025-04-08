function [SkyMag,magnitude] = Extinction(SkyMag_ori, magnitude_ori, zenith, wave)
    % This function is used to calculate the conversion between the magnitude 
    % corresponding to the actual number of star photons received and the 
    % actual magnitude based on the atmospheric extinction ratio.
    
    wavelist = [0.36,0.43,0.55,0.65,0.82];
    k_list = [0.55,0.25,0.15,0.09,0.06];
    temme = abs(wave-wavelist);
    [~,nwav] = min(temme);
    k = k_list(nwav);
    SkyMag = SkyMag_ori + secd(zenith)*k;
    magnitude = magnitude_ori + secd(zenith)*k;
end