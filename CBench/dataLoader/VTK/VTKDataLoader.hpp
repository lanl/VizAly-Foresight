/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#ifndef _VTK_LOADER_H_
#define _VTK_LOADER_H_

#include <sstream>
#include <string>
#include <unordered_map>

#include <vtkPointData.h>
#include <vtkMPIController.h>
#include <vtkTrivialProducer.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>

#include <vtkXMLPImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>

#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "dataLoaderInterface.hpp"

#include "json.hpp"
#include "timer.hpp"
#include "utils.hpp"


class VTKDataLoader: public DataLoaderInterface
{
	vtkSmartPointer<vtkImageData> imageData;
	vtkSmartPointer<vtkImageData> imageWriteData;

	vtkSmartPointer<vtkMPIController> contr;

	int numRanks;
	int myRank;	

	std::string elementType; // cell or point
	int numComponents;

  public:
	VTKDataLoader();
	~VTKDataLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

    int saveInputFileParameters() { return 1; };
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value){};
  	bool loadUncompressedFields(nlohmann::json const&) { return false; } 
};

inline VTKDataLoader::VTKDataLoader()
{
	imageData = nullptr;
	
	contr = vtkSmartPointer<vtkMPIController>::New();
	contr->Initialize(NULL, NULL, 1);
}


inline VTKDataLoader::~VTKDataLoader()
{

}


inline void VTKDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);

	// Read in file and store in imageData
	vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  	reader->SetFileName(filename.c_str());
  	reader->Update();

  	imageData = reader->GetOutput();
	
  	debugLog << "vtkXMLImageDataReader Number of point array: " << reader->GetNumberOfPointArrays() << std::endl;
  	debugLog << "vtkXMLImageDataReader Number of cell array: " << reader->GetNumberOfCellArrays() << std::endl;
}



inline int VTKDataLoader::loadData(std::string paramName)
{
	Timer clock("load");

	debugLog << "vtkImageData Number of points: " << imageData->GetNumberOfPoints() << std::endl;
    debugLog << "vtkImageData Number of cells: " << imageData->GetNumberOfCells() << std::endl;
	debugLog << "Number of cell array: " << imageData->GetCellData()->GetNumberOfArrays() << std::endl;
	debugLog << "Number of point array: " << imageData->GetPointData()->GetNumberOfArrays() << std::endl;


	// Get Dimensions
	int *dims = imageData->GetDimensions();
	

	// Read in array in myArray
	vtkSmartPointer<vtkDataArray> myArray;

	bool found = false;
	for (int i=0; i<imageData->GetCellData()->GetNumberOfArrays(); i++)
	{
		myArray = imageData->GetCellData()->GetArray( i );
		std::string arrName(myArray->GetName());
		if (arrName == paramName)
		{
			elementType = "cell";
			numComponents = myArray->GetNumberOfComponents();

			// Cells is dims -1
			origDims[0] = dims[0]-1;
			origDims[1] = dims[1]-1;
			origDims[2] = dims[2]-1;

			// Get type
			elemSize = vtkDataArray::GetDataTypeSize(myArray->GetDataType());
			switch (myArray->GetDataType())
			{
				case VTK_FLOAT:
					dataType = "float";
					break;
				case VTK_DOUBLE:
					dataType = "double";
					break;
				case VTK_INT:
					dataType = "int";
					break;
				default:
					if (myRank == 0)
						std::cout << "This datatype is not supported yet for VTK-VTK. This program will now exit!" << std::endl;

					MPI_Finalize();
			}

			found = true;
			break;
		}
	}


	if (!found)
	{
		for (int i=0; i<imageData->GetPointData()->GetNumberOfArrays(); i++)
		{
			myArray = imageData->GetPointData()->GetArray( i );
			std::string arrName(myArray->GetName());
			if (arrName == paramName)
			{
				elementType = "point";
				numComponents = myArray->GetNumberOfComponents();

				origDims[0] = dims[0];
				origDims[1] = dims[1];
				origDims[2] = dims[2];

				// Get type
				elemSize = vtkDataArray::GetDataTypeSize(myArray->GetDataType());
				switch (myArray->GetDataType())
				{
					case VTK_FLOAT:
						dataType = "float";
						break;
					case VTK_DOUBLE:
						dataType = "double";
						break;
					case VTK_INT:
						dataType = "int";
						break;
				}

				std::cout << "Found point array " << paramName <<  std::endl;
				break;
			}
		}
	}


	if (!found)
	{
		std::cout << "This field does not exist!!!" << std::endl;
		return 0;
	}



	totalNumberOfElements = origDims[0]*origDims[1]*origDims[2];	

	//
	// TODO: MPI split!!!
	Partition current = getPartition(myRank, numRanks, origDims[0], origDims[1], origDims[2]);
	rankOffset[0] = current.min_x;
	rankOffset[1] = current.min_y;
	rankOffset[2] = current.min_z;

	sizePerDim[0] = current.max_x - current.min_x;
	sizePerDim[1] = current.max_y - current.min_y;
	sizePerDim[2] = current.max_z - current.min_z;

	numElements = sizePerDim[0]*sizePerDim[1]*sizePerDim[2];



	// Create space for data and store data there
	allocateMem(dataType, numElements, 0, data);


	size_t indexLocal = 0;
	for (size_t z=0; z<sizePerDim[2]; z++)
		for (size_t y=0; y<sizePerDim[1]; y++)
			for (size_t x=0; x<sizePerDim[0]; x++)
			{
				size_t indexGlobal = (rankOffset[2] + z) * (origDims[0] * origDims[1]) +
									 (rankOffset[1] + y) *  origDims[0] +
									 (rankOffset[0] + x);

				// Choose based on datatype
				if (dataType == "float")
					((float *)data)[indexLocal] = myArray->GetComponent(indexGlobal, 0);
				else if (dataType == "double")
					((double *)data)[indexLocal] = myArray->GetComponent(indexGlobal, 0);
				else if (dataType == "int")
					((int *)data)[indexLocal] = myArray->GetComponent(indexGlobal, 0);
				indexLocal++;
			}

	

	// Sould we save the data
	if (saveData)
	{
		imageWriteData = vtkSmartPointer<vtkImageData>::New();

		double spacing[3];
		imageData->GetSpacing (spacing);
		imageWriteData->SetSpacing (spacing[0], spacing[1], spacing[2]);

		double origin[3];
		imageData->GetOrigin(origin);
		imageWriteData->SetOrigin (origin[0], origin[1], origin[2]);
	}
	

	clock.stop("load");
	debugLog << "origDims: "   << origDims[0]   << ", " << origDims[1]   << ", " << origDims[2]   << std::endl;
	debugLog << "sizePerDim: " << sizePerDim[0] << ", " << sizePerDim[1] << ", " << sizePerDim[2] << std::endl;
	debugLog << "rankOffset: " << rankOffset[0] << ", " << rankOffset[1] << ", " << rankOffset[2] << std::endl;
	debugLog << "numElements: " << numElements << std::endl;
	debugLog << "totalNumberOfElements: " << totalNumberOfElements << std::endl;
	
	debugLog << "Loading data took: " << clock.getDuration("load") << " s" << std::endl;
}



inline int VTKDataLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock("save");

		

	if (dataType == "float")
	{
		vtkSmartPointer<vtkFloatArray> temp = vtkFloatArray::New();

		temp->SetNumberOfTuples(numElements);
		temp->SetNumberOfComponents(numComponents);
		temp->SetName(paramName.c_str());


		for (size_t index=0; index<numElements; index++)
			temp->SetValue( index, ((float *)cData)[index]);

		imageWriteData->GetCellData()->AddArray(temp);
	}
	else if (dataType == "double")
	{
		vtkSmartPointer<vtkDoubleArray> temp = vtkDoubleArray::New();

		temp->SetNumberOfTuples(numElements);
		temp->SetNumberOfComponents(numComponents);
		temp->SetName(paramName.c_str());


		for (size_t index=0; index<numElements; index++)
			temp->SetValue( index, ((double *)cData)[index]);

		imageWriteData->GetCellData()->AddArray(temp);
	}
	else if (dataType == "int")
	{
		vtkSmartPointer<vtkIntArray> temp = vtkIntArray::New();

		temp->SetNumberOfTuples(numElements);
		temp->SetNumberOfComponents(numComponents);
		temp->SetName(paramName.c_str());


		for (size_t index=0; index<numElements; index++)
			temp->SetValue( index, ((int *)cData)[index]);

		imageWriteData->GetCellData()->AddArray(temp);
	}
	




	imageWriteData->SetDimensions(sizePerDim[0], sizePerDim[1], sizePerDim[2]);
	imageWriteData->SetExtent(rankOffset[0], rankOffset[0]+sizePerDim[0],
							  rankOffset[1], rankOffset[1]+sizePerDim[1],
							  rankOffset[2], rankOffset[2]+sizePerDim[2]);


	clock.stop("save");
	debugLog << "saving data took: " << clock.getDuration("save") << " s" << std::endl;
}



inline int VTKDataLoader::writeData(std::string _filename)
{
	Timer clock("write");


	vtkSmartPointer<vtkXMLPImageDataWriter> writer = vtkSmartPointer<vtkXMLPImageDataWriter>::New();
	std::string outputFilename;

  	if (numRanks > 1)
	{
		_filename.replace((_filename.length()-3),3,"pvti");  
		outputFilename =   _filename;

		vtkNew<vtkTrivialProducer> tp;
		tp->SetOutput(imageWriteData);
		tp->SetWholeExtent(0, origDims[0],
		                   0, origDims[1],
		                   0, origDims[2]);

		writer->SetInputConnection(tp->GetOutputPort());
		writer->SetController(contr);
	}
	else
		outputFilename = _filename + ".vti";



	writer->SetDataModeToBinary();
    writer->SetCompressor(nullptr);
    writer->SetWriteSummaryFile(1);
	writer->SetFileName(outputFilename.c_str());
	writer->SetNumberOfPieces(numRanks);
	writer->SetStartPiece(myRank);
    writer->SetEndPiece(myRank);


	#if VTK_MAJOR_VERSION <= 5
    	writer->SetInput(imageWriteData);
	#else
  		writer->SetInputData(imageWriteData);
	#endif

	writer->Write();
  	

  	clock.stop("write");
	debugLog << "writing data took: " << clock.getDuration("write") << " s" << std::endl;
}


#endif

