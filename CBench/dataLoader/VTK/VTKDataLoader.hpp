/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Jesus Pulido
================================================================================*/

#ifndef _VTK_LOADER_H_
#define _VTK_LOADER_H_

#include <sstream>
#include <string>
#include <unordered_map>

#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkFloatArray.h>

#include "dataLoaderInterface.hpp"

#include "json.hpp"
#include "timer.hpp"
#include "utils.hpp"


class VTKDataLoader: public DataLoaderInterface
{
	vtkSmartPointer<vtkImageData> imageData;
	vtkSmartPointer<vtkImageData> imageWriteData;

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
	


	

  	log << "vtkXMLImageDataReader Number of point array: " << reader->GetNumberOfPointArrays() << std::endl;
  	log << "vtkXMLImageDataReader Number of cell array: " << reader->GetNumberOfCellArrays() << std::endl;
}



inline int VTKDataLoader::loadData(std::string paramName)
{
	Timer clock;
	clock.start();

	log << "vtkImageData Number of points: " << imageData->GetNumberOfPoints() << std::endl;
    log << "vtkImageData Number of cells: " << imageData->GetNumberOfCells() << std::endl;
	log << "Number of cell array: " << imageData->GetCellData()->GetNumberOfArrays() << std::endl;
	log << "Number of point array: " << imageData->GetPointData()->GetNumberOfArrays() << std::endl;


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
		exit(0);
	}


	totalNumberOfElements = origDims[0]*origDims[1]*origDims[2];	
	sizePerDim[0] = origDims[0];
	sizePerDim[1] = origDims[1];
	sizePerDim[2] = origDims[2];


	//
	// TODO: MPI split!!!
	Partition current = getPartition(myRank, numRanks, origDims[0], origDims[1], origDims[2]);
	numElements = totalNumberOfElements;
	//


	// Create space for data and store data there
	allocateMem("float", numElements, 0, data);

	//data = myArray->GetVoidPointer(0);
	for (int i=0; i<numElements; i++)
		((float *)data)[i] =  myArray->GetComponent(i, 0);
	


	imageWriteData = vtkSmartPointer<vtkImageData>::New();
	imageWriteData->SetDimensions(origDims[0]+1, origDims[1]+1, origDims[2]+1);

	
	clock.stop();
	log << "Loading data took: " << clock.getDuration() << " s" << std::endl;
}



inline int VTKDataLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock;
	clock.start();


	// TODO: MPI Gather to a new array
	void *allData;
	if (myRank == 0)
		allocateMem(dataType, totalNumberOfElements, 0, allData, );
		

	
		vtkSmartPointer<vtkFloatArray> temp = vtkFloatArray::New();;

		temp->SetNumberOfTuples(numElements);
		temp->SetNumberOfComponents(numComponents);
		temp->SetName(paramName.c_str());

		temp->SetNumberOfValues(numElements);
		for (size_t i=0; i<numElements; i++)
			temp->SetValue( i, ((float *)cData)[i]);

		imageWriteData->GetCellData()->AddArray(temp);
	
	if (myRank == 0)
		deAllocateMem(dataType, allData);

	clock.stop();
	log << "saving data took: " << clock.getDuration() << " s" << std::endl;
}



inline int VTKDataLoader::writeData(std::string _filename)
{
	Timer clock;
	clock.start();

	if (myRank == 0)
	{
		vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	  	writer->SetFileName(_filename.c_str());
		#if VTK_MAJOR_VERSION <= 5
	  		writer->SetInputConnection(imageWriteData->GetProducerPort());
		#else
	  		writer->SetInputData(imageWriteData);
		#endif
	  	writer->Write();
  	}

  	clock.stop();
	log << "writing data took: " << clock.getDuration() << " s" << std::endl;
}


#endif
