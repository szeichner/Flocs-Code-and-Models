classdef Photo
    %PHOTO This class converts photos to grayscale and populates different
    %matrix properties that define values associated with the picture
    
    properties
        PixelMatrixAvgVector
        pixelMatrix
    end
    
    methods
        function obj = Photo(photoDirectory, fileName, photoCropArray)
            %PHOTO Construct an instance of this class
            %   create a photo object from a file name, and populate a
            %   matrix and vector that describes average pixel values
            %   across rows of the photo matrix
            
            %read file, convert to gray scale and simplify data format for
            %matrix
            photoImport = imread(char(strcat(photoDirectory, '/', fileName)));
            
            photoImport = mat2gray(photoImport);
            
            %truncate photo to get just strip of interest
            photoImport = imcrop(photoImport, photoCropArray);

            obj.pixelMatrix = single(squeeze(photoImport(:,:,1)));
            [~, columns] = size(obj.pixelMatrix());
         
            %take average value of each row of matrix and put into vector
            PixelMatrixAvgVector = zeros(columns, 1);
            for i= 1:size(obj.pixelMatrix,1)
                PixelMatrixAvgVector(i) = mean(obj.pixelMatrix(i,:));
            end
            
            obj.PixelMatrixAvgVector = PixelMatrixAvgVector;
            
        end
        
    end
end

