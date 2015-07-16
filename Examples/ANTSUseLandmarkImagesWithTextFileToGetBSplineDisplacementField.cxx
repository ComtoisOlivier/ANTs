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
	if (file.is_open())
	{
		while (file.good())
		{
			getline(file, line);
			char * charLine = new char[line.length() + 1];
			strcpy(charLine, line.c_str());
			token = std::strtok(charLine, ",");
			while (token != NULL)
			{
				result.push_back(std::atof(token));
				token = strtok(NULL, ",");
			}
			delete[] charLine;
		}
		file.close();
	}
	else
	{
		throw std::runtime_error(std::string("Cannot open file : ")+filename);
	}
	return result;
}

template<unsigned int ImageDimension>
void GetRealValuePointSetFromFile(typename itk::PointSet<float, ImageDimension>::Pointer &curved,typename itk::PointSet<float, ImageDimension>::Pointer &straight,string curvedFilename, string straightFilename) {
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
	typedef float RealType;
	typedef unsigned int LabelType;
	typedef itk::Image<LabelType, ImageDimension> LabelImageType;
	typedef itk::Image<RealType, ImageDimension> ReadImageType;

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

	typedef itk::ImportImageFilter<LabelType, ImageDimension> ImporterType;
	typename ImporterType::Pointer importer = ImporterType::New();
	importer->SetImportPointer(	const_cast<LabelType *>(fixedImage->GetBufferPointer()), numberOfPixels, filterHandlesMemory);
	importer->SetRegion(fixedImage->GetBufferedRegion());
	importer->SetOrigin(fixedImage->GetOrigin());
	importer->SetSpacing(fixedImage->GetSpacing());
	importer->SetDirection(identityDirection);
	importer->Update();
	const typename ImporterType::OutputImageType * parametricInputImage = importer->GetOutput();

	typename ImageReaderType::Pointer movingReader = ImageReaderType::New();
	movingReader->SetFileName(argv[2]);
	movingReader->Update();
	typename LabelImageType::Pointer movingImage = movingReader->GetOutput();

///////////////////////////////////////////////////////////////////////////////////////////

	typedef itk::Vector<RealType, ImageDimension> VectorType;
	typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

	typedef itk::PointSet<LabelType, ImageDimension> PointSetType;
	typename PointSetType::Pointer fixedPoints = PointSetType::New();
	fixedPoints->Initialize();

	std::vector<LabelType> fixedLabels;

	itk::ImageRegionIteratorWithIndex<LabelImageType> ItF(fixedImage,
			fixedImage->GetLargestPossibleRegion());
//
	typename PointSetType::Pointer movingPoints = PointSetType::New();
	movingPoints->Initialize();
//
	std::vector<LabelType> movingLabels;

	itk::ImageRegionIteratorWithIndex<LabelImageType> ItM(movingImage,
			movingImage->GetLargestPossibleRegion());


	std::vector<float> labelWeights;
	std::vector<LabelType> userLabels;

	//initialize points container
	typedef itk::PointSet<RealType, ImageDimension> RealPointSetType;
	typename RealPointSetType::Pointer fixedPts = RealPointSetType::New();
	typename RealPointSetType::Pointer movingPts = RealPointSetType::New();

	fixedPts->Initialize();
	movingPts->Initialize();

	//Checking if both optional inputs are present
	if (string(argv[6]).compare("") != 0) {
		if (string(argv[7]).compare("") != 0) {
			try {
				string straightFilename = argv[7];
				string curvedFilename = argv[6];
				GetRealValuePointSetFromFile<3>(movingPts, fixedPts,
						curvedFilename, straightFilename);
			} catch (exception &e) {
				throw;
			}
		} else
			throw std::runtime_error(
					"Both arguments the straight and curved file must be specified");
	} else
		throw std::runtime_error(
				"A file must be specified in order to take input real labels.");
	// Now match up the center points

	typedef itk::PointSet<VectorType, ImageDimension> DisplacementFieldPointSetType;
	typedef itk::BSplineScatteredDataPointSetToImageFilter<DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
	typedef typename BSplineFilterType::WeightsContainerType WeightsContainerType;

	typename WeightsContainerType::Pointer weights = WeightsContainerType::New();
	weights->Initialize();
	const typename WeightsContainerType::Element weight = 1.0;

	typename DisplacementFieldPointSetType::Pointer fieldPoints = DisplacementFieldPointSetType::New();
	fieldPoints->Initialize();
	unsigned long count = 0;

	typename RealPointSetType::PointsContainerConstIterator mIt = movingPts->GetPoints()->Begin();
//	typename PointSetType::PointsContainerConstIterator mItI = movingCenters->GetPoints()->Begin();
	typename RealPointSetType::PointDataContainerIterator mItD = movingPts->GetPointData()->Begin();
	typename RealPointSetType::PointsContainerConstIterator fIt = fixedPts->GetPoints()->Begin();
//	typename PointSetType::PointsContainerConstIterator fItI = fixedCenters->GetPoints()->Begin();
	typename RealPointSetType::PointDataContainerIterator fItD = fixedPts->GetPointData()->Begin();

	//transform points to physcial index
	typedef itk::PointSet<RealType, ImageDimension> RealPhysicalPoints;
//	unsigned int pointId = 0;

	typename RealPhysicalPoints::Pointer realMovingPts = RealPhysicalPoints::New();
	typename RealPhysicalPoints::Pointer realFixedPts = RealPhysicalPoints::New();
	realFixedPts->Initialize();
	realMovingPts->Initialize();

//	}
	while (mItD != movingPts->GetPointData()->End()) {
		fIt = fixedPts->GetPoints()->Begin();
		fItD = fixedPts->GetPointData()->Begin();
		while (fItD != fixedPts->GetPointData()->End()) {
			if (fItD.Value() == mItD.Value()) {
				typename RealPointSetType::PointType fpoint = fIt.Value();
				typename RealPointSetType::PointType mpoint = mIt.Value();

				VectorType vector;

				typename LabelImageType::PointType indexPoint;
				typename DisplacementFieldType::PointType fixedPhysicalPoint;
				typename DisplacementFieldType::PointType movingPhysicalPoint;
				typename DisplacementFieldType::PointType fieldPoint;
				itk::ContinuousIndex<double, ImageDimension> fixedCidx;
				itk::ContinuousIndex<double, ImageDimension> movingCidx;
				for (unsigned int i = 0; i < ImageDimension; i++) {
					fixedCidx[i] = fpoint[i];
					movingCidx[i] = mpoint[i];
					//vector[i] = mpoint[i] - fpoint[i];
					//fixedCidx[i] = indexPoint[i];
				}

				fixedImage->TransformContinuousIndexToPhysicalPoint(fixedCidx, fixedPhysicalPoint);
				movingImage->TransformContinuousIndexToPhysicalPoint(movingCidx, movingPhysicalPoint);

				for (unsigned int i = 0; i < ImageDimension; i++) {

					vector[i] = movingPhysicalPoint[i] - fixedPhysicalPoint[i];

				}

				fixedImage->TransformPhysicalPointToContinuousIndex( fixedPhysicalPoint, fixedCidx );

				parametricInputImage->TransformContinuousIndexToPhysicalPoint(fixedCidx, fieldPoint);

				printf("\nFIELD POINT : x = %f   y = %f   z = %f", fieldPoint[0], fieldPoint[1], fieldPoint[2]);

				fieldPoints->SetPoint(count, fieldPoint);
				fieldPoints->SetPointData(count, vector);

				weights->InsertElement(count, weight);
				count++;

				break;
			}
			++fItD;
			++fIt;
		}

		++mItD;
		++mIt;
	}

	bool enforceStationaryBoundary = true;
	if (argc > 9) {
		enforceStationaryBoundary = static_cast<bool>(atoi(argv[9]));
	}
	if (enforceStationaryBoundary) {
		typename LabelImageType::IndexType startIndex2 =
				fixedImage->GetLargestPossibleRegion().GetIndex();

		typename LabelImageType::SizeType inputSize2 = fixedImage->GetLargestPossibleRegion().GetSize();
		for (ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF) {
			typename LabelImageType::IndexType index = ItF.GetIndex();

			bool isOnStationaryBoundary = false;
			for (unsigned int d = 0; d < ImageDimension; d++) {
				if (index[d] == startIndex2[d] || index[d] == startIndex2[d] + static_cast<int>(inputSize2[d]) - 1)
				{
					isOnStationaryBoundary = true;
					break;
				}
			}

			if (isOnStationaryBoundary) {
				VectorType vector;

				vector.Fill(0.0);

				typename PointSetType::PointType fixedPoint;
				parametricInputImage->TransformIndexToPhysicalPoint(index,
						fixedPoint);

				fieldPoints->SetPoint(count, fixedPoint);
				fieldPoints->SetPointData(count, vector);
				count++;
			}
		}
	}
	typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

	unsigned int numberOfLevels = atoi(argv[5]);

	unsigned int splineOrder = 3;
	if (argc > 8) {
		splineOrder = atoi(argv[8]);
	}

	std::vector<unsigned int> meshSize = ConvertVector<unsigned int>(std::string(argv[4]));
	typename BSplineFilterType::ArrayType ncps;
	ncps.Fill(0);


	if (meshSize.size() == 1) {
		ncps.Fill(meshSize[0] + splineOrder);
	} else if (meshSize.size() == ImageDimension) {
		for (unsigned int d = 0; d < ImageDimension; d++) {
			ncps[d] = meshSize[d] + splineOrder;
		}
	} else {
		std::cerr << "Invalid meshSize format." << std::endl;
	}
//   std::cout << ncps << std::endl;
//
//    bspliner->DebugOn();
	bspliner->SetOrigin(fixedImage->GetOrigin());
	bspliner->SetSpacing(fixedImage->GetSpacing());
	bspliner->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	bspliner->SetDirection(fixedImage->GetDirection());
	bspliner->SetGenerateOutputImage(true);
	bspliner->SetNumberOfLevels(numberOfLevels);
	bspliner->SetSplineOrder(splineOrder);
	bspliner->SetNumberOfControlPoints(ncps);
	bspliner->SetInput(fieldPoints);
	bspliner->SetPointWeights(weights);
	bspliner->Update();

	typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType,
			RealType> InterpolatorType;
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetInputImage(bspliner->GetOutput());

	std::cout << std::endl << "Distance errors:" << std::endl;

	mIt = movingPts->GetPoints()->Begin();
	mItD = movingPts->GetPointData()->Begin();

	while (mItD != movingPts->GetPointData()->End())
	{
	    fIt = fixedPts->GetPoints()->Begin();
	    fItD = fixedPts->GetPointData()->Begin();

		while (fItD != fixedPts->GetPointData()->End())
		{
			if (fItD.Value() == mItD.Value())
			{
				itk::ContinuousIndex<double, ImageDimension> fpointIdx;
				itk::ContinuousIndex<double, ImageDimension> mpointIdx;
				typename DisplacementFieldType::PointType fpoint;
				typename DisplacementFieldType::PointType mpoint;
				typename RealPointSetType::PointType nonPhysicalFPoint = fIt.Value();
				typename RealPointSetType::PointType nonPhysicalMPoint = mIt.Value();

				for (unsigned int i = 0; i < ImageDimension; i++)
				{
					fpointIdx[i] = nonPhysicalFPoint[i];
					mpointIdx[i] = nonPhysicalMPoint[i];
				}

				fixedImage->TransformContinuousIndexToPhysicalPoint(fpointIdx, fpoint);
				movingImage->TransformContinuousIndexToPhysicalPoint(mpointIdx, mpoint);

				printf("\nMPOINT : x = %f   y = %f   z = %f\n", mpoint[0], mpoint[1], mpoint[2]);
				printf("\nFPOINT : x = %f   y = %f   z = %f\n", fpoint[0], fpoint[1], fpoint[2]);

				VectorType displacement = (mpoint - fpoint);
				printf("\nDisplacement : %f, %f, %f\n", displacement[0], displacement[1], displacement[2]);

				typename InterpolatorType::PointType point;
				for (unsigned int i = 0; i < ImageDimension; i++) {
					point[i] = fpoint[i];
				}
				printf("\nEvaluating point :: x = %f   y = %f   z = %f\n", point[0], point[1], point[2]);
				//VectorType vector = interpolator->Evaluate(point);
				VectorType vector;
				if(interpolator->IsInsideBuffer(point))
				{
					vector = interpolator->Evaluate(point);
				}
				else
				{
					throw std::runtime_error("Point is not inside image bounds");
				}
				printf("\nCALCULATING ERROR\n");
				RealType error = (vector - displacement).GetNorm();
				std::cout << "  " << fItD.Value() << ": " << error << std::endl;

				break;
			}
			++fItD;
			++fIt;
		}

		++mItD;
		++mIt;
	}

	typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(argv[3]);
	writer->SetInput(bspliner->GetOutput());
	writer->Update();
	return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ANTSUseLandmarkImagesWithTextFileToGetBSplineDisplacementField(
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
