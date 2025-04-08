 classdef Observer
    %It is used to define the entire observer of Telescope array
    %   The exposure calculator for telescope arrays
    %   It holds all telescopes in a telescope array and also some
    %   important parameters for exposure calculations.
    
    properties
        Rearth_m =  6371.23e3; % earth radius in metre
        atmh_m = 50e3; %Total atmospheric height in metre
        iPhotometryzero = [0.36 397500;0.43 1077120;0.55 860000;0.65 953610;...
            0.82 659400];
        %iWaveDataMatrix is the wavelength matrix along with effectiveness
        %Other Zero point reference
        % Photometric system from Rigaut's OOMAO
        %  Band ( wavelength , bandwidth , zero point )
        %         U  ( 0.360e-6 , 0.070e-6 , 2.0e12 )
        %         B  ( 0.440e-6 , 0.100e-6 , 5.4e12 )
        %         V0  ( 0.500e-6 , 0.090e-6 , 3.3e12 )
        %         V  ( 0.550e-6 , 0.090e-6 , 3.3e12 )
        %         R  ( 0.640e-6 , 0.150e-6 , 4.0e12 )
        %         I  ( 0.790e-6 , 0.150e-6 , 2.7e12 )
        %         J  ( 1.215e-6 , 0.260e-6 , 1.9e12 )
        %         H  ( 1.654e-6 , 0.290e-6 , 1.1e12 )
        %         Ks ( 2.157e-6 , 0.320e-6 , 5.5e11 )
        %         K  ( 2.179e-6 , 0.410e-6 , 7.0e11 )
        %         L  ( 3.547e-6 , 0.570e-6 , 2.5e11 )
        %         M  ( 4.769e-6 , 0.450e-6 , 8.4e10 )
        %         Na ( 0.589e-6 , 0        , 3.3e12 )
        %         EOS ( 1.064e-6 , 0        , 3.3e12 )
        %         HK  ( (1.654e-6*1.1e12+2.179e-6*7.0e11)/(1.1e12+7.0e11)
        %           , 0.290e-6+00.410e-6 , 1.1e12+7.0e11 )
        %SDSS 0.367 75065;0.4686 546600;0.6166 493485; 0.7480 359665
        Efocal = 5; %Focal length of the telescope in metre
        Aperture = 1.0; %Aperture size of the telescope in metre
        Fieldofview = 3.0; %Field of view defintion in degree
        TelEfficiency = 0.76; %Total efficiency of the telescope ref AST3
        Camsize = 8*1e-6; %Size of pixel of the CCD in mm
        CamPixelNumber = 2048; %number of pixels in CCD
        CamDarkCurrent = 1; %Dark current in e-/s
        CamReadOutNoise = 1; %Readout noise in e-
        CamWellDepth = 1e7; %Full well depth for simulation in e-
        CamType = 'CCD'; %Type of the camera
        CamPrice = 1e8; %It is used to define the price of camera
        CamEff = 0.99; %Quantum Efficiency of the camera
        SeidelCoeff  = [0, 0.6361,0.0063, 0.7061, -0.9986, 0.8118];

        Fov = 1;%Default Field of view in degree
        Wavelength = 0.6231;%Default wavelength for simulation in um
        SkybackNoise  = 22.3; 

        iWaveDataMatrix = [0.43, 1; 0.55 1; 0.65, 1];
        %The fov distribution is defined by a matrix
        fovdistribution = {};
        %Ntel is the number of telescopes used in the array
        Ntel = 75;
    end
    
    methods
        function Observer = Observer(Efocal, Aperture, TelEfficiency, CamSize,...
                CamPixelNumber,CamDarkCurrent,CamReadOutNoise,CamWellDepth,...
                SeidelCoeff,Fov,Wavelength,SkybackNoise)
            if nargin > 0
                Observer.Efocal = Efocal; 
                Observer.Aperture = Aperture; 
                Observer.TelEfficiency = TelEfficiency; 
                Observer.CamSize = CamSize;
                Observer.CamPixelNumber = CamPixelNumber; 
                Observer.CamDarkCurrent=CamDarkCurrent; 
                Observer.CamReadOutNoise = CamReadOutNoise;
                Observer.CamWellDepth = CamWellDepth; 
                Observer.SkybackNoise = SkybackNoise; 
                Observer.SeidelCoeff  = SeidelCoeff; 
                Observer.Fov = Fov;
                Observer.Wavelength = Wavelength; 
            end
        end
    end
    
    methods
        %%%%%Function that are used to define a telescope that would be
        %%%%%used in Sitian
        function [psfcube,Observer] = definetelescope(Observer,TheApplication, args,Aperture,Fov)
            %Define numbers of pixels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 修改 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Get the pixelscale of camera in arcsec
            temcampixelscale = Observer.Camsize/Observer.Efocal*206265;
            %We could obtain the fov of this camera
            temfov = temcampixelscale*1.414*Observer.CamPixelNumber/3600;
            maxfov = Observer.Fov/2;
            rdiv = 2;
            thetadiv = 6;
            fovusage = [maxfov,rdiv,thetadiv];
            [psfcube,Observer] = Observer.deploystate(TheApplication, args,fovusage);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [simulatePSF] = simulateImageByPSF(Observer, psfdata)
            %Get the pixelscale of camera in arcsec
            temcampixelscale = Observer.Camsize/Observer.Efocal*206265;
        end
        %%%%%%%Realistic PSFs are generated from PSF cube and atmospheric
        %%%%%%%turbulence
        function psfcube = realisticpsf(Observer,psfcube,seeingfwhm)
            %We generate realistic psfcube from this part
            Npsf = size(psfcube{1});
            focallength = Observer.Efocal;
            for ind = 1:Npsf(2)
                psfcube{1}{ind} = Observer.onerealisticpsf(...
                    psfcube{1}{ind},psfcube{2}(ind),psfcube{3}(ind),...
                    focallength,seeingfwhm);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%One realistic PSF is used to generate realistic PSF
        function finalpsf = onerealisticpsf(Observer,onepsf,pixelscale,...
                wavelength,focallength,seeingfwhm)
            %We will set the seeingfwhm according to wavelength
            scale = (wavelength/0.606)^(-0.2);
            sigma = seeingfwhm/2.3548*scale;
            %We get the size of psf
            pixelscalearcsec = pixelscale*1e-6/focallength*206265;

            psf = onepsf;
            newsize = size(psf);

            atmpsf = fspecial('gaussian',newsize(1),sigma/pixelscalearcsec);
            finalpsf = conv2(psf,atmpsf,'same');
            finalpsf = finalpsf/sum(finalpsf(:));
            finalpsf(isnan(finalpsf)) = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Observer] = rescal_lens(Observer, TheApplication, args, Aperture, workingfolder, telescopefolder, telescopename, Fov)
            import ZOSAPI.*;
            TheSystem = TheApplication.PrimarySystem;
            sampleDir = TheApplication.SamplesDir;
            TheAnalyses = TheSystem.Analyses;
            % Add your custom code here...--------------------
            path1 = pwd;
            [pathstr,name,suffix] = fileparts(telescopename);
            workingfolder = strcat(path1,'\',workingfolder);

            if ~exist(workingfolder, 'dir') 
                mkdir(workingfolder); 
            else
%                 disp('exist'); 
            end

            outfile = strcat(path1,'\workingfolder\',name,'.zmx');
            if contains(pathstr, ':') 
                orgfile = strcat(telescopename);
%                 orgfile = telescopename;
                fprintf('working file path：%s \n',orgfile);
            else
                orgfile = strcat(path1, '\telescopetemplate\', name, '.zmx');
            end

            try
                result = TheSystem.LoadFile(orgfile, false);
                TheSystem.SaveAs(outfile);
            catch
                fprintf('copy error!');
            end

            SystemApertureData = TheSystem.SystemData.Aperture.ApertureValue;
            ApertureValue = SystemApertureData/1000;   
            Value = Aperture/ApertureValue;

            NewField1 = TheSystem.SystemData.Fields.AddField(0, Fov/2, 1.00);
            sysField = TheSystem.SystemData.Fields;
            sysField.SetFieldType(ZOSAPI.SystemData.FieldType.Angle);

            ScaleLens = TheSystem.Tools.OpenScale(); 
            ScaleLens.ScaleByFactor = true;
            ScaleLens.ScaleFactor = Value;
            ScaleLens.RunAndWaitForCompletion();
            ScaleLens.Close();

            Observer.Aperture = Aperture;
            TheMFE = TheSystem.MFE;
            Operand_1 = TheMFE.GetOperandAt(1);
            Operand_1.ChangeType(ZOSAPI.Editors.MFE.MeritOperandType.EFFL);
            TheMFE.CalculateMeritFunction();
            EFFL = Operand_1.Value;

            Observer.Efocal = EFFL*10^-3;

            TheSystem.Save();
            % -------------------------------------------
%             r = [];
        end

        
        function [psfcube,Observer] = deploystate(Observer, TheApplication, args,fovusage)
            %In this code, should have wavelength set.

            maxfov = fovusage(1);  
            rdiv = fovusage(2); 
            thetadiv = fovusage(3); 
            fovusagetemr = linspace(0,maxfov,rdiv+6);   
            fovusagetemtheta = linspace(0,360,thetadiv-4);  
            Nr = size(fovusagetemr,2);  
            Ndir = size(fovusagetemtheta,2); 
            fovusagetemtheta = fovusagetemtheta(1:end-1);
            psfmatrix = {};
            pixelscale=[];
            wavelist=[];
            Nwav = size(Observer.iWaveDataMatrix);
            Nwav = Nwav(1);
            temfovdistribution = {};
            for indw = 1:Nwav
                for indr =1:Nr    %2
                    for indtheta = 1:Ndir-1    % 6
                        temx=fovusagetemr(indr)*cosd(fovusagetemtheta(indtheta));
                        temy=fovusagetemr(indr)*sind(fovusagetemtheta(indtheta));
                        wavelength = Observer.iWaveDataMatrix(indw,1);
                        temdata = Observer.getpsf(TheApplication, args,temx,temy,wavelength);
                        psfmatrix{end+1} = temdata{1};
                        pixelscale = [pixelscale,temdata{2}];
                        wavelist = [wavelist,wavelength];
                        temfovdistribution{end+1} = {temx,temy,wavelength};
                        if fovusagetemr(indr) == 0
                            break
                        end
                    end
                end
            end
            Observer.fovdistribution = temfovdistribution;
            psfcube = {psfmatrix,pixelscale,wavelist};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%Function to get psf for a particular field of view and
        %%%%%%wavelength
        function psfme = getpsf(Observer,TheApplication, args, fovx,fovy,wavelength)
            %We need to set the field of view parameters
            %We can set the field of view
            import ZOSAPI.*;
            TheSystem = TheApplication.PrimarySystem;
            TheAnalyses = TheSystem.Analyses;
            % Add your custom code here...--------------------
            sysField = TheSystem.SystemData.Fields;
            sysField.GetField(1).X = fovx;
            sysField.GetField(1).Y = fovy;
            sysField.GetField(1).Weight = 1;
            %As well we set the wave
            sysWave = TheSystem.SystemData.Wavelengths;
            sysWave.GetWavelength(1).Wavelength = wavelength;
            % set analyses setting
            newWin = TheAnalyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.HuygensPsf);
            newWin_setting = newWin.GetSettings();
            newWin_setting.Wavelength.SetWavelengthNumber(1);    
            newWin_setting.Field.SetFieldNumber(1);
            newWin_setting.Type = ZOSAPI.Analysis.Settings.HuygensPsfTypes.Linear;
            newWin_setting.UseCentroid = true;  
            newWin_setting.PupilSampleSize = ZOSAPI.Analysis.SampleSizes.S_32x32;
            newWin_setting.ImageSampleSize = ZOSAPI.Analysis.SampleSizes.S_64x64;
            newWin.ApplyAndWaitForCompletion();
            newWin_Results = newWin.GetResults();
            filepath = strcat(pwd, '\PSFAnalysis1.txt');
            
            newWin_Results.GetTextFile(filepath);
            TheAnalyses.CloseAnalysis(newWin);
            TheSystem.Save();

            psf_data = importdata(filepath);
            psf = psf_data.data;
            temstr = psf_data.textdata{6};
            numreg = regexp(temstr,'-?\d*\.?\d*','match');
            pixelscale = str2double(numreg{1});
            psfme= {psf,pixelscale};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%Function used to obtain final images for a camera
        function Observer = definecam(Observer,Camera)
            %Observer is the class definition
            %There are 6 cameras here:
            %Type PixelNumber PixelSize Price        Dark   RON   Eff  Fullw
            %CMOS1 8k          10um     1860000      0.0011  1    0.95  50000
            %CMOS2 4k          9um       120000      0.0011  1    0.95  50000
            %CCD1  4k          15um     1300000      0.001   1   0.97  160000
            %CCD2  6k          15um     2200000      0.001   1   0.97  160000
            %CCD3  9k          10um     3550000      0.001   1   0.97  160000
            %CCD4  10k          9um     4000000      0.001   1   0.97  160000
            if strcmp(Camera,'CCD')
                Observer.Camsize = 9*1e-6;
                Observer.CamEff = 0.97;
                Observer.CamDarkCurrent = 0.0011;
                Observer.CamReadOutNoise = 1;
                Observer.CamWellDepth = 160000;
                Observer.CamType = 'CCD';
            elseif strcmp(Camera,'CMOS')
                Observer.Camsize = 9*1e-6;
                Observer.CamEff = 0.95;
                Observer.CamDarkCurrent = 0.0011;
                Observer.CamReadOutNoise = 1;
                Observer.CamWellDepth = 50000;
                Observer.CamType = 'CMOS';
            else
                disp('No such device!')
            end
        end

        function [imagecube,photons] = getimagecube(Observer,psfcube,mag,wavelength,...
                pixelsize,ron,dark,exposuretime,efficiency,stackerror,D_radio)
            %This function would consider noise, response and other things
            campixel = pixelsize;
            %We obtain the wavelength
            %Pick the observation wavelength band
            temme = abs(wavelength-Observer.iWaveDataMatrix);
            Nwav = size(Observer.iWaveDataMatrix);
            Nwav = Nwav(1);
            [~,Index] = min(temme(:,1));
            wavelength = Observer.iWaveDataMatrix(Index,:);
            %We load the psfcube with their perspective pixelscale
            TotalNpsfcube = size(psfcube{1});
            TotalNpsfcube = TotalNpsfcube(2);
            %We will process these PSFs to generate real observation data
            psfmat = psfcube{1};
            psfpixelscale = psfcube{2};
            %Image cube is the cube for all data for different fov and
            %different wavelength
            %Arranged as the same size of psfcube
            imagecube={};
            photons={};
            Npsf = TotalNpsfcube/Nwav;
            for ind = 1:Npsf
                temid = (Index-1)*Npsf+ind;
                tempsf = psfmat{temid};
                psf_size = size(tempsf);  
                tempixelscale = psfpixelscale(temid)*1e-6; 

                if stackerror > 1e-10 && exposuretime > 30*60
                    temstackerror = stackerror/(tempixelscale/Observer.Efocal*...
                        206265);
                    tempsf = imgaussfilt(tempsf,temstackerror); 

                    tempsf(tempsf<0) = 0;
                    tempsf(isnan(tempsf)) = 0;
                end
                scalefactor = tempixelscale/campixel;
                tempsf = imresize(tempsf,scalefactor);
                tempsf(isnan(tempsf)) = 0;
                tempsf(tempsf<0)=0;
                tempsf = tempsf/sum(tempsf(:));
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

                imagecube{end+1} = {outimg,noisemat,orgimg};
                photons{end+1} = {allphotons,skyback};
            end
            function [outimg,noisemat,orgimg,allphotons,skyback]=photonnum(exposuretime)
                sizeOfPSF = size(tempsf);
                allphotons = magzero*exposuretime*pi*...
                    (Observer.Aperture*100/2)^2*D_radio*10^(-1*mag/2.5)*efficiency;
                orgimg = poissrnd(allphotons*(tempsf/sum(sum(tempsf(:)))));
                %Calculate the sky background noise
                skyback = magzero*exposuretime*0.5...
                     *(Observer.Aperture*100)^2*D_radio*10^(-1*Observer.SkybackNoise/2.5)*efficiency*0.5*((Observer.Fov*3600)^2);
                skyback = skyback/((Observer.CamPixelNumber)^2);
                noisemat = poissrnd(skyback,size(tempsf))+...
                    poissrnd(dark*exposuretime,size(tempsf))+...
                    poissrnd((ron^2),size(tempsf));   
                outimg = noisemat+orgimg;
                outimg(outimg>Observer.CamWellDepth) ...   
                    = Observer.CamWellDepth;
                noisemat(noisemat>Observer.CamWellDepth) =...
                    Observer.CamWellDepth;
                orgimg(orgimg>Observer.CamWellDepth) =...
                    Observer.CamWellDepth;
            end
        end

        % Function to evaluate SNR
        function [snrmat,posmat,distance_half_,signalphoton_,noisephoton_,signalphotons,skybackphotons] = ...
                calsnrmap(Observer,imgcube, psfcubere, wavelength,camerapixcelsize,photons)
            %It is a function used to calcualte distribution of SNR in the
            Nimg = size(Observer.fovdistribution);
            Nimg = Nimg(2);
            Nwav = size(Observer.iWaveDataMatrix);
            Nwav = Nwav(1);
            temme = abs(wavelength-Observer.iWaveDataMatrix);
            [~,nwav] = min(temme(:,1));
            totalimg = Nimg/Nwav;   
            posmat = ones(2,totalimg);
            snrmat = ones(1,totalimg);
            distance_half_ = ones(1,totalimg);
            signalphoton_ = ones(1,totalimg);
            noisephoton_ = ones(1,totalimg);
            signalphotons = ones(1,totalimg);
            skybackphotons = ones(1,totalimg);
            for ind = 1:totalimg
                temind = (nwav-1)*totalimg+ind;

                sigmat = imgcube{ind}{3};
                noisemat = imgcube{ind}{2};
                signalphotons(ind) = photons{ind}{1};
                skybackphotons(ind) = photons{ind}{2};

                psf_afterconv = psfcubere{1}{temind};
                pixelscale_afterconv = psfcubere{2}(temind);
                [snrmat(ind),distance_half_(ind),signalphoton_(ind),noisephoton_(ind)] = ...
                                calsnr11(sigmat,noisemat,psf_afterconv,pixelscale_afterconv,camerapixcelsize);
         
                posmat(1,ind) = Observer.fovdistribution{temind}{1};
                posmat(2,ind) = Observer.fovdistribution{temind}{2};
            end
        end

        function totalprice = CalTotalPrice(Observer,Ntel,Aperture,exposuretime)
            %We first evalute the price for telescope build
            telprice = Observer.CalTelPrice(Aperture);
            Maintainceprice = telprice*0.15;
            camprice = Observer.CamPrice;
            totalprice=(telprice+Maintainceprice+camprice)*Ntel;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Function used to calculate the price of a telescope
        function telprice = CalTelPrice(Observer, Aperture)
            %Based on two price:
            %The Price is Proportional to the aperture
            %The Price is 12000000 for a 1 metre with 7 degree
            %The Price is 6000000 for a 1 metre with 3 degree
            fovlist = [3,7];
            fovprice = [6000000,12000000];
            onemprice = interp1(fovlist,fovprice,Observer.Fov,'linear','extrap');
            %Scale to parameter
            telprice = onemprice*(Aperture)^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [NUM,allfov,radio_Field] = CalNumOfEnough(Observer,Overlap,overallfov,width_height)
            % Used to calculate how many telescopes are needed to complete the target area in the current field of view.
            fovonetel = 0.5*(Observer.Fov)^2;
            NUM = ceil(overallfov/fovonetel);
            k = (NUM/width_height(2)/width_height(1))^0.5;
            m = width_height(1)*k;   
            n = width_height(2)*k;    

            radio_Field = (2*2^0.5)*Observer.Fov/(Observer.CamPixelNumber*(Observer.Camsize/Observer.Efocal/pi*180));
            wastedline = Overlap*(Observer.Camsize/Observer.Efocal/pi*180)*radio_Field;

            wastedfov = ((m*(n-1))+(n*(m-1)))*wastedline*(2*2^0.5)*Observer.Fov-3*(m-1)*(n-1)*wastedline^2;            
            allfov_ = 0.5*(Observer.Fov)^2*NUM-wastedfov;

            while allfov_ < overallfov
                NUM = NUM+1;
                wastedfov = ((m*(n-1))+(n*(m-1)))*wastedline*(2*2^0.5)*Observer.Fov-3*(m-1)*(n-1)*wastedline^2; 
                allfov_ = 0.5*(Observer.Fov)^2*NUM-wastedfov;
            end
             allfov = allfov_;
        end
    end
end