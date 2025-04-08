function [CamPrice, pixcel_num_width] = CalCamPriceZOS(TheApplication, args, camtype, pixelsize, price, height, width)
    % This function calculates the approximate price of the corresponding CCD or CMOS camera based on the size of the image plane in the optical system.
    % camtype            input: Camera type, can be either CCD or CMOS, default type is CCD ------ CCD:0, CMOS:1
    % pixelsize          input: Pixel size of the corresponding camera sensor (in mm), default pixel size is 9 microns, 0.009 mm
    % price              input: Camera price per square millimeter, default price is the CCD price of 494/mm^2
    % height             input: Aspect ratio corresponding to the height of the camera sensor, default aspect ratio is 1
    % width              input: Aspect ratio corresponding to the width of the camera sensor
    
    % CamPrice           output: Corresponding camera price
    % pixcel_num_width   output: Number of pixels in width
    
    import ZOSAPI.*;
    TheSystem = TheApplication.PrimarySystem;
    sampleDir = TheApplication.SamplesDir;
    TheAnalyses = TheSystem.Analyses;

    if nargin == 7
        K = height/width;
    elseif nargin == 5
        K = 1;
    elseif nargin == 4
        K = 1;
        price = caprice(camtype);
    elseif nargin == 3
        pixelsize = 9/1000;
        price = caprice(camtype);
        K = 1;
    elseif nargin == 2
        pixelsize = 9/1000;
        price = 494;
        K = 1;
    else
        error('please ener the correct number of parameters')
    end

    try
        NUM = TheSystem.LDE.NumberOfSurfaces;
    catch ErrorInfo
        disp(ErrorInfo);
        disp('没有正常获取系统面数')
    end

    Surface = TheSystem.LDE.GetSurfaceAt(NUM-1);
    SemiDiaOfSurf = Surface.SemiDiameter;

    width_num = floor((((SemiDiaOfSurf^2)/(K^2 +1))^(0.5))*2/pixelsize); 
    height_num = floor((K * width_num));                                
    pixcel_num_width = width_num;                                    
    width_ = pixelsize * width_num;                                    
    height_ = pixelsize * height_num;                                
    CamPrice = width_ * height_ *  price;                                

    function price = caprice(type)   
        % type = 0： CCD， 1： CMOS
        if type == 1
            price = 494;
        elseif type == 2
            price = 247;
        else
            error('please ener the correct parameters')
        end
    end

end
     
