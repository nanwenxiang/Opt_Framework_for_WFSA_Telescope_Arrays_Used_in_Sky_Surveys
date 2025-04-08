function [imagecube,photons] = getimagecube_replace_single(Observer,psfdata,samplingscale,mag,wavelength,...
    pixelsize,ron,dark,exposuretime,efficiency,stackerror,D_radio)
    % Used for processing PSF images, simulating imaging through photonization.

    campixel = pixelsize;
    %We obtain the wavelength
    temme = abs(wavelength-Observer.iWaveDataMatrix);
    [~,Index] = min(temme(:,1));
    % We will process these PSFs to generate real observation data
    psfpixelscale = samplingscale;
    imagecube={};
    photons={};
    psf_size = size(psfdata); 
    tempixelscale = psfpixelscale*1e-6;
    if stackerror > 1e-10 && exposuretime > 30*60
        temstackerror = stackerror/(tempixelscale/Observer.Efocal*...
            206265);
        psfdata = imgaussfilt(psfdata,temstackerror); 

        psfdata(psfdata<0) = 0;
        psfdata(isnan(psfdata)) = 0;
    end
    scalefactor = tempixelscale/campixel;
    psfdata = imresize(psfdata,scalefactor);
    psfdata(isnan(psfdata)) = 0;
    psfdata(psfdata<0)=0;
    psfdata = psfdata/sum(psfdata(:));
    magzero = Observer.iPhotometryzero(Index,2);
    
    if exposuretime < 30*60
        [outimg,noisemat,orgimg,allphotons,skyback] = photonnum(exposuretime);
    else  
        outimg=0;noisemat=0;orgimg=0;allphotons=0;skyback=0;
        for num = 1:floor(exposuretime/1800)
            [outimg_,noisemat_,orgimg_,allphotons_,skyback_] = photonnum(30*60);
            outimg=outimg_+outimg;noisemat=noisemat_+noisemat;orgimg=orgimg_+orgimg;
            allphotons=allphotons+allphotons_;skyback=skyback+skyback_;
        end
        [outimg_,noisemat_,orgimg_,allphotons_,skyback_] = photonnum(rem(exposuretime,1800)); 
        outimg=outimg_+outimg;noisemat=noisemat_+noisemat;orgimg=orgimg_+orgimg;
        allphotons=allphotons+allphotons_;skyback=skyback+skyback_;
    end

    outimg = poissrnd(outimg,size(outimg));
    imagecube{end+1} = {outimg,noisemat,orgimg};
    photons{end+1} = {allphotons,skyback};


    function [outimg,noisemat,orgimg,allphotons,skyback]=photonnum(exposuretime)
        sizeOfPSF = size(psfdata);
        allphotons = magzero*exposuretime*pi*...
            (Observer.Aperture*100/2)^2*D_radio*10^(-1*mag/2.5)*efficiency;
        orgimg = allphotons*(psfdata/sum(sum(psfdata(:))));
        %Calculate the sky background noise
        skyback = magzero*exposuretime*0.5...
             *(Observer.Aperture*100)^2*D_radio*10^(-1*Observer.SkybackNoise/2.5)*efficiency*0.5*((Observer.Fov*3600)^2);
        skyback = skyback/((Observer.CamPixelNumber)^2);
        noisemat = skyback+...
                dark*exposuretime+...
                (ron^2);  
        outimg = noisemat+orgimg;
        sigma = 25; 
        gaussian_noise = sigma * randn(size(outimg)); 
        outimg_with_noise = outimg + gaussian_noise;
        outimg(outimg>Observer.CamWellDepth) ...   
            = Observer.CamWellDepth;
        noisemat(noisemat>Observer.CamWellDepth) =...
            Observer.CamWellDepth;
        orgimg(orgimg>Observer.CamWellDepth) =...
            Observer.CamWellDepth;
    end
end