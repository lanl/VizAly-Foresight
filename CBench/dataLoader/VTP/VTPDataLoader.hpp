/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#ifndef _VTP_LOADER_H_
#define _VTP_LOADER_H_

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

#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>


#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "dataLoaderInterface.hpp"

#include "json.hpp"
#include "timer.hpp"
#include "utils.hpp"


class VTPDataLoader: public DataLoaderInterface
{
	vtkSmartPointer<vtkPolyData> polyData;
	vtkSmartPointer<vtkPolyData> polyWriteData;

	vtkSmartPointer<vtkMPIController> contr;

	int numRanks;
	int myRank;	

	std::string elementType; // cell or point
	int numComponents;

  public:
	VTPDataLoader();
	~VTPDataLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

    int saveInputFileParameters() { return 1; };
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value){};
  	bool loadUncompressedFields(nlohmann::json const&) { return false; } 
};

inline VTPDataLoader::VTPDataLoader()
{
	polyData = nullptr;
	
	contr = vtkSmartPointer<vtkMPIController>::New();
	contr->Initialize(NULL, NULL, 1);
}


inline VTPDataLoader::~VTPDataLoader()
{

}


inline void VTPDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);

	// Read in file and store in imageData
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  	reader->SetFileName(filename.c_str());
  	reader->Update();

  	polyData = reader->GetOutput();
	
  	debugLog << "vtkXMLPolyDataReader Number of point array: " << reader->GetNumberOfPointArrays() << std::endl;
  	debugLog << "vtkXMLPolyDataReader Number of cell array: " << reader->GetNumberOfCellArrays() << std::endl;
}



inline int VTPDataLoader::loadData(std::string paramName)
{
	Timer clock("load");

	debugLog << "vtkPolyData Number of points: " << polyData->GetNumberOfPoints() << std::endl;
    debugLog << "vtkPolyData Number of cells: " << polyData->GetNumberOfCells() << std::endl;
	debugLog << "Number of cell array: " << polyData->GetCellData()->GetNumberOfArrays() << std::endl;
	debugLog << "Number of point array: " << polyData->GetPointData()->GetNumberOfArrays() << std::endl;


	// Get Dimensions (physical range)
	//int *dims = polyData->GetDimensions();
	

	// Read in array in myArray
	vtkSmartPointer<vtkDataArray> myArray;

	bool found = false;
	//CellData
	for (int i=0; i<polyData->GetCellData()->GetNumberOfArrays(); i++)
	{
		myArray = polyData->GetCellData()->GetArray( i );
		std::string arrName(polyData->GetCellData()->GetArrayName(i));
		if (arrName == paramName)
		{
			elementType = "cell";
			numComponents = myArray->GetNumberOfComponents();

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

	// PointData
	if (!found)
	{
		for (int i=0; i<polyData->GetPointData()->GetNumberOfArrays(); i++)
		{
			myArray = polyData->GetPointData()->GetArray( i );
			std::string arrName(polyData->GetPointData()->GetArrayName(i));
			if (arrName == paramName)
			{
				elementType = "point";
				numComponents = myArray->GetNumberOfComponents();

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
				found = true;
				break;
			}
		}
	}

    // If xyz, convert to I
	// Convert to I and use lossless? lossy + quant?

	if (!found)
	{
		std::cout << "This field does not exist!!!" << std::endl;
		return 0;
	}



	//totalNumberOfElements = origDims[0]*origDims[1]*origDims[2];	
	totalNumberOfElements = polyData->GetNumberOfPoints();
	std::cout << "Points: " << totalNumberOfElements << std::endl;

	//
	// TODO: MPI split!!!
	Partition current = getPartition(myRank, numRanks, origDims[0], origDims[1], origDims[2]);
	rankOffset[0] = current.min_x;
	rankOffset[1] = current.min_y;
	rankOffset[2] = current.min_z;

	sizePerDim[0] = current.max_x - current.min_x;
	sizePerDim[1] = current.max_y - current.min_y;
	sizePerDim[2] = current.max_z - current.min_z;

	//numElements = sizePerDim[0]*sizePerDim[1]*sizePerDim[2];

	numElements = totalNumberOfElements;


	// Create space for data and store data there
	allocateMem(dataType, numElements, 0, data);

	sizePerDim[0] = numElements;	// For compression

	/*size_t indexLocal = 0;
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

	*/
	for (size_t i=0; i<numElements; i++)
	{
		if (dataType == "float")
			((float *)data)[i] = myArray->GetComponent(i, 0);
		else if (dataType == "double")
			((double *)data)[i] = myArray->GetComponent(i, 0);
		else if (dataType == "int")
			((int *)data)[i] = myArray->GetComponent(i, 0);
	}

	// Sould we save the data
	if (saveData)
	{
		polyWriteData = vtkSmartPointer<vtkPolyData>::New();
		// Pass through x/y/z positions
		inOutDataPoints = polyData->GetPoints();

		double spacing[3];
		//polyData->GetSpacing (spacing);
		//polyWriteData->SetSpacing (spacing[0], spacing[1], spacing[2]);

		double origin[3];
		//polyData->GetOrigin(origin);
		//polyWriteData->SetOrigin (origin[0], origin[1], origin[2]);
	}
	

	clock.stop("load");
	debugLog << "origDims: "   << origDims[0]   << ", " << origDims[1]   << ", " << origDims[2]   << std::endl;
	debugLog << "sizePerDim: " << sizePerDim[0] << ", " << sizePerDim[1] << ", " << sizePerDim[2] << std::endl;
	debugLog << "rankOffset: " << rankOffset[0] << ", " << rankOffset[1] << ", " << rankOffset[2] << std::endl;
	debugLog << "numElements: " << numElements << std::endl;
	debugLog << "totalNumberOfElements: " << totalNumberOfElements << std::endl;
	
	debugLog << "Loading data took: " << clock.getDuration("load") << " s" << std::endl;
	return 1; // All good!
}



inline int VTPDataLoader::saveCompData(std::string paramName, void * cData)
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

		polyWriteData->GetPointData()->AddArray(temp);
	}
	else if (dataType == "double")
	{
		vtkSmartPointer<vtkDoubleArray> temp = vtkDoubleArray::New();

		temp->SetNumberOfTuples(numElements);
		temp->SetNumberOfComponents(numComponents);
		temp->SetName(paramName.c_str());


		for (size_t index=0; index<numElements; index++)
			temp->SetValue( index, ((double *)cData)[index]);

		polyWriteData->GetPointData()->AddArray(temp);
	}
	else if (dataType == "int")
	{
		vtkSmartPointer<vtkIntArray> temp = vtkIntArray::New();

		temp->SetNumberOfTuples(numElements);
		temp->SetNumberOfComponents(numComponents);
		temp->SetName(paramName.c_str());


		for (size_t index=0; index<numElements; index++)
			temp->SetValue( index, ((int *)cData)[index]);

		polyWriteData->GetPointData()->AddArray(temp);
	}
	polyWriteData->SetPoints(inOutDataPoints);

	//polyWriteData->SetDimensions(sizePerDim[0], sizePerDim[1], sizePerDim[2]);
	//polyWriteData->SetExtent(rankOffset[0], rankOffset[0]+sizePerDim[0],
	//						  rankOffset[1], rankOffset[1]+sizePerDim[1],
	//						  rankOffset[2], rankOffset[2]+sizePerDim[2]);


	clock.stop("save");
	debugLog << "saving data took: " << clock.getDuration("save") << " s" << std::endl;
	
}



inline int VTPDataLoader::writeData(std::string _filename)
{
	Timer clock("write");


	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	std::string outputFilename;

	outputFilename = _filename + ".vtp";



	writer->SetDataModeToBinary();
    writer->SetCompressor(nullptr);
    //writer->SetWriteSummaryFile(1);
	writer->SetFileName(outputFilename.c_str());
	//writer->SetNumberOfPieces(numRanks);
	//writer->SetStartPiece(myRank);
    //writer->SetEndPiece(myRank);


	#if VTK_MAJOR_VERSION <= 5
    	writer->SetInput(polyWriteData);
	#else
  		writer->SetInputData(polyWriteData);
	#endif

	writer->Write();
  	

  	clock.stop("write");
	debugLog << "writing data took: " << clock.getDuration("write") << " s" << std::endl;
}


#endif

