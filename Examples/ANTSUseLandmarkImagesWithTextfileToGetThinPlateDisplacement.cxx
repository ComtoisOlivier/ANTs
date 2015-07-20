/** ANTS Landmarks used to initialize an b-spline displacement field ... */

#include "antsUtilities.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImportImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkPointSet.h"
#include <itkImageFunction.h>
#include "itkThinPlateSplineKernelTransform.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>
#include <string.h>
#include <vector>
#include <stdlib.h>

using namespace std;

namespace ants {

std::vector<float> ReadLabelValueFromFile(std::string filename) {
	std::ifstream file(filename.c_str());
	string line;
	char* token;
	std::vector<float> result;
	std::vector<float>::iterator it;

	it = result.begin();
	if (file.is_open()) {
		while (file.good()) {
			getline(file, line);
			char * charLine = new char[line.length() + 1];
			strcpy(charLine, line.c_str());
			token = std::strtok(charLine, ",");
			while (token != NULL) {
				result.push_back(std::atof(token));
				token = strtok(NULL, ",");
			}
			delete[] charLine;
		}
		file.close();
	} else {
		throw std::runtime_error(std::string("Cannot open file : ") + filename);
	}
	return result;
}

template<unsigned int ImageDimension>
void GetRealValuePointSetFromFile(
		typename itk::PointSet<float, ImageDimension>::Pointer &curved,
		typename itk::PointSet<float, ImageDimension>::Pointer &straight,
		string curvedFilename, string straightFilename) {
	//define variables needed for this function
	string line;
	string token;
	std::vector<float> straightValues;
	std::vector<float> curvedValues;

	//fetching coordinates from files
	try {
		straightValues = ReadLabelValueFromFile(straightFilename);
		curvedValues = ReadLabelValueFromFile(curvedFilename);
		//check if not empty
		if (curvedValues.size() == 0) {
			throw std::runtime_error(
					"No value found in the curved labels text file.");
		}
		if (straightValues.size() == 0) {
			throw std::runtime_error(
					"No value found in straight the labels text file.");
		}
		if ((straightValues.size() % 3) != 0) {
			throw std::runtime_error(
					"The straight file size must be divisible by 3 (since there is 3 dimensions to the image)");
		}
		if ((curvedValues.size() % 3) != 0) {
			throw std::runtime_error(
					"The curved file size must be divisible by 3 (since there is 3 dimensions to the image)");
		}
		if (straightValues.size() != curvedValues.size()) {
			throw std::runtime_error(
					"Curved and Straight files do not have the same size.");
		}
	} catch (const std::exception &e) {
		throw;
	}

	//Produce the Straight points
	typedef float Realtype;
	typedef itk::PointSet<Realtype, ImageDimension> PointsFromText;

	std::vector<float>::iterator itS = straightValues.begin();
	std::vector<float>::iterator itC = curvedValues.begin();
	unsigned int pointId = 0;

	while (itS != straightValues.end()) {
		typename PointsFromText::PointType pC;
		typename PointsFromText::PointType pS;
		for (int i = 0; i < 3; i++) {
			pS[i] = *itS;
			pC[i] = *itC;
			itS++;
			itC++;
		}
		//Place points in the point sets
		curved->SetPoint(pointId, pC);
		curved->SetPointData(pointId, pointId);
		straight->SetPoint(pointId, pS);
		straight->SetPointData(pointId, pointId++);
	}
	printf("RETURNING POINT VALUES");
	return;

}

template<unsigned int ImageDimension>
int LandmarkBasedWithTextDisplacementFieldTransformInitializer(int argc, char *argv[]) {

	////////////////////////////////////////////////////////////////////////////////////////////
	/*CREATING USEFUL TYPES*/
	////////////////////////////////////////////////////////////////////////////////////////////
	typedef float RealType;
	typedef unsigned int LabelType;
	typedef itk::Image<LabelType, ImageDimension> LabelImageType;
	typedef itk::Image<RealType, ImageDimension> ReadImageType;
	typedef itk::Vector< float, ImageDimension > FieldVectorType;
	typedef itk::PointSet<LabelType, ImageDimension> PointSetType;
	typedef itk::ResampleImageFilter<LabelImageType, LabelImageType> ResamplerType;
	typedef itk::LinearInterpolateImageFunction<LabelImageType, double > InterpolatorType;
	typedef itk::ImageFileWriter< LabelImageType >  DeformedImageWriterType;
	typedef itk::Image< FieldVectorType,  ImageDimension > DisplacementFieldType;
	typedef itk::ImageFileWriter< DisplacementFieldType > FieldWriterType;
	////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////
	/*LOADING FIXED IMAGE*/
	////////////////////////////////////////////////////////////////////////////////////////////
	typedef itk::ImageFileReader<LabelImageType> ImageReaderType;
	typename ImageReaderType::Pointer fixedReader = ImageReaderType::New();
	fixedReader->SetFileName(argv[1]);
	fixedReader->Update();
	typename LabelImageType::Pointer fixedImage = fixedReader->GetOutput();
	typename LabelImageType::DirectionType fixedDirection =	fixedImage->GetDirection();
	typename LabelImageType::DirectionType fixedDirectionInverse(fixedDirection.GetInverse());
	typename LabelImageType::DirectionType identityDirection; //No clue
	identityDirection.SetIdentity();

	const typename LabelImageType::RegionType & bufferedRegion = fixedImage->GetBufferedRegion();
	const itk::SizeValueType numberOfPixels = bufferedRegion.GetNumberOfPixels();

	const bool filterHandlesMemory = false;
	////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////
	/*CREATING IMPORTER FOR FIXED IMAGE*/
	////////////////////////////////////////////////////////////////////////////////////////////
	typedef itk::ImportImageFilter<LabelType, ImageDimension> ImporterType;
	typename ImporterType::Pointer importer = ImporterType::New();
	importer->SetImportPointer(const_cast<LabelType *>(fixedImage->GetBufferPointer()),	numberOfPixels, filterHandlesMemory);
	importer->SetRegion(fixedImage->GetBufferedRegion());
	importer->SetOrigin(fixedImage->GetOrigin());
	importer->SetSpacing(fixedImage->GetSpacing());
	importer->SetDirection(identityDirection);
	importer->Update();
	const typename ImporterType::OutputImageType * parametricInputImage = importer->GetOutput();
	////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////
	/*LOADING MOVING IMAGE*/
	////////////////////////////////////////////////////////////////////////////////////////////
	typename ImageReaderType::Pointer movingReader = ImageReaderType::New();
	movingReader->SetFileName(argv[2]);
	movingReader->Update();
	typename LabelImageType::Pointer movingImage = movingReader->GetOutput();

	typedef itk::Vector<RealType, ImageDimension> VectorType;
	typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
	////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////
	/*INIT IMAGE CONTAINERS*/
	////////////////////////////////////////////////////////////////////////////////////////////
	itk::ImageRegionIteratorWithIndex<LabelImageType> ItF(fixedImage, fixedImage->GetLargestPossibleRegion());
	itk::ImageRegionIteratorWithIndex<LabelImageType> ItM(movingImage, movingImage->GetLargestPossibleRegion());

	//initialize points container
	typedef itk::PointSet<RealType, ImageDimension> RealPointSetType;
	typename RealPointSetType::Pointer fixedPts = RealPointSetType::New();
	typename RealPointSetType::Pointer movingPts = RealPointSetType::New();

	fixedPts->Initialize();
	movingPts->Initialize();
	////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////
	/*IMPORTING LANDMARKS FROM TEXT FILES*/
	////////////////////////////////////////////////////////////////////////////////////////////
	if (string(argv[6]).compare("") != 0) {
		if (string(argv[7]).compare("") != 0) {
			try {
				string straightFilename = argv[7];
				string curvedFilename = argv[6];
				GetRealValuePointSetFromFile<3>(movingPts, fixedPts, curvedFilename, straightFilename);
			} catch (exception &e) {
				throw;
			}
		} else
			throw std::runtime_error(
					"Both arguments the straight and curved file must be specified");
	} else
		throw std::runtime_error(
				"A file must be specified in order to take input real labels.");
	////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////
	/*DEFINING KERNEL TRANSFORM*/
	////////////////////////////////////////////////////////////////////////////////////////////
	typedef itk::Vector< float, ImageDimension > FieldVectorType;
	typedef itk::Image< FieldVectorType,  ImageDimension > DisplacementFieldType;
	typedef itk::ThinPlateSplineKernelTransform< RealType, ImageDimension> TransformType;
	TransformType::Pointer tps = TransformType::New();

	tps->SetSourceLandmarks(movingPts);
	tps->SetTargetLandmarks(fixedPts);
	tps->ComputeWMatrix();
	////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////////////////////
	/*DEFINING RESAMPLER TO EXECUTE TRANSFORM*/
	////////////////////////////////////////////////////////////////////////////////////////////
	LabelImageType::ConstPointer inputImage = fixedReader->GetOutput();
	ResamplerType::Pointer resampler = ResamplerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	resampler->SetInterpolator( interpolator );

    InputImageType::SpacingType spacing = inputImage->GetSpacing();
    InputImageType::PointType   origin  = inputImage->GetOrigin();
    InputImageType::DirectionType direction  = inputImage->GetDirection();
    InputImageType::RegionType region = inputImage->GetBufferedRegion();
    InputImageType::SizeType   size =  region.GetSize();

    resampler->SetOutputSpacing( spacing );
    resampler->SetOutputDirection( direction );
    resampler->SetOutputOrigin(  origin  );
    resampler->SetSize( size );
    resampler->SetTransform( tps );

    resampler->SetOutputStartIndex(  region.GetIndex() );
    resampler->SetInput( fixedReader->GetOutput() );

    DeformedImageWriterType::Pointer deformedImageWriter = DeformedImageWriterType::New();
    deformedImageWriter->SetInput( resampler->GetOutput() );
    deformedImageWriter->SetFileName( argv[3] );

    try
    {
       deformedImageWriter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
       std::cerr << "Exception thrown " << std::endl;
       std::cerr << excp << std::endl;
       return EXIT_FAILURE;
    }

    DisplacementFieldType::Pointer field = DisplacementFieldType::New();
    field->SetRegions( region );
    field->SetOrigin( origin );
    field->SetSpacing( spacing );
    field->Allocate();

    typedef itk::ImageRegionIterator< DisplacementFieldType > FieldIterator;
    FieldIterator fi( field, region );
    fi.GoToBegin();
    TransformType::InputPointType  point1;
    TransformType::OutputPointType point2;
    DisplacementFieldType::IndexType index;
    itk::ContinuousIndex<double, ImageDimension> fixedCidx;
    itk::ContinuousIndex<double, ImageDimension> movingCidx;

    FieldVectorType displacement;
    while( ! fi.IsAtEnd() )
    {
      index = fi.GetIndex();
      field->TransformIndexToPhysicalPoint( index, point1 );
      point2 = tps->TransformPoint( point1 );
      for ( unsigned int i = 0;i < ImageDimension;i++)
      {
        displacement[i] = point2[i] - point1[i];
      }
      fi.Set( displacement );
      ++fi;
      }
    //Write computed deformation field
    FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName( argv[4] );
    fieldWriter->SetInput( field );
    try
      {
      fieldWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
	////////////////////////////////////////////////////////////////////////////////////////////

}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ANTSUseLandmarkImagesWithTextFileToGetThinPlateDisplacementField(
		std::vector<std::string> args, std::ostream* /*out_stream = NULL */) {
	// put the arguments coming in as 'args' into standard (argc,argv) format;
	// 'args' doesn't have the command name as first, argument, so add it manually;
	// 'args' may have adjacent arguments concatenated into one argument,
	// which the parser should handle
	args.insert(args.begin(),
			"ANTSUseLandmarkImagesWithTextFileToGetBSplineDisplacementField");
	int argc = args.size();
	char* * argv = new char *[args.size() + 1];
	for (unsigned int i = 0; i < args.size(); ++i) {
		// allocate space for the string plus a null character
		argv[i] = new char[args[i].length() + 1];
		std::strncpy(argv[i], args[i].c_str(), args[i].length());
		// place the null character in the end
		argv[i][args[i].length()] = '\0';
	}
	argv[argc] = ITK_NULLPTR;
	// class to automatically cleanup argv upon destruction
	class Cleanup_argv {
	public:
		Cleanup_argv(char* * argv_, int argc_plus_one_) :
				argv(argv_), argc_plus_one(argc_plus_one_) {
		}

		~Cleanup_argv() {
			for (unsigned int i = 0; i < argc_plus_one; ++i) {
				delete[] argv[i];
			}
			delete[] argv;
		}

	private:
		char* * argv;
		unsigned int argc_plus_one;
	};
	Cleanup_argv cleanup_argv(argv, argc + 1);

	// antscout->set_stream( out_stream );

	if (argc < 4) {
		std::cerr << "Usage:   " << argv[0]
				<< " fixedImageWithLabeledLandmarks  movingImageWithLabeledLandmarks outputDisplacementField "
				<< "meshSize[0]xmeshSize[1]x... numberOfLevels movingLabelsFilename fixedLabelsFilename [order=3] [enforceStationaryBoundaries=1]"
				<< std::endl;
		std::cerr
				<< " we expect the input images to be (1) N-ary  (2) in the same physical space as the images you want to "
				<< std::endl;
		std::cerr
				<< " register and (3 ) to have the same landmark points defined within them ... "
				<< std::endl;
		std::cerr
				<< " landmarks will be defined from the center of mass of the labels in the input images . "
				<< std::endl;
		std::cerr << " You can use ITK-snap to generate the label images. "
				<< std::endl;
		std::cerr
				<< " The optional landmarks weights are read from a text file where each row is either:"
				<< std::endl;
		std::cerr << " \"label,labelWeight\" or " << std::endl;
		std::cerr << " \"labelWeight\" or " << std::endl;
		std::cerr
				<< " If the latter format is used, the label weights are assumed to be arranged in ascending order by label."
				<< std::endl;
		if (argc >= 2
				&& (std::string(argv[1]) == std::string("--help")
						|| std::string(argv[1]) == std::string("-h"))) {
			return EXIT_SUCCESS;
		}
		return EXIT_FAILURE;
	}

	// Get the image dimension
	std::string fn = std::string(argv[1]);
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
			fn.c_str(), itk::ImageIOFactory::ReadMode);
	imageIO->SetFileName(fn.c_str());
	imageIO->ReadImageInformation();

	switch (imageIO->GetNumberOfDimensions()) {
	case 2: {
		//LandmarkBasedWithTextDisplacementFieldTransformInitializer<2>(argc, argv);
		throw std::runtime_error("This function only takes 3D images");
	}
		break;
	case 3: {
		try {
			LandmarkBasedWithTextDisplacementFieldTransformInitializer<3>(argc,
					argv);
		} catch (const std::exception &e) {
			std::cerr << e.what() << std::endl;
		}
	}
		break;
	default:
		std::cerr << "Unsupported dimension" << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
} // namespace ants
