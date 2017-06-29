/*
 * Copyright © 2017 Arno Mayrhofer
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * http://www.wtfpl.net/ for more details.
 *
 */

#include <stdio.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPImageDataWriter.h>
#include <vtkImageData.h>
#include <mpi.h>

#include <vtkImageAlgorithm.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkDataArray.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkAlgorithmOutput.h>
#include <vtkMPIController.h>

class ImageSource : public vtkImageAlgorithm
{
public:
    ImageSource();
    void setWholeExtent(const int * const extent);
    void setOrigin(const double * const origin);
    void setSpacing(const double * const spacing);

private:
    double origin_[3], spacing_[3];
    int wholeExtent_[6];

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *outputVector);

    void ExecuteDataWithInformation(vtkDataObject *, vtkInformation *outInfo);

    int FillOutputPortInformation(int port, vtkInformation *info);
};

ImageSource::ImageSource()
{
    origin_[0] = 0.0;
    origin_[1] = 0.0;
    origin_[2] = 0.0;
    spacing_[0] = 1.0;
    spacing_[1] = 1.0;
    spacing_[2] = 1.0;
    wholeExtent_[0] = 0;
    wholeExtent_[1] = 1;
    wholeExtent_[2] = 0;
    wholeExtent_[3] = 1;
    wholeExtent_[4] = 0;
    wholeExtent_[5] = 1;

    this->SetNumberOfInputPorts(0);
}

void ImageSource::setWholeExtent(const int * const extent)
{
    for (int i = 0; i < 6; i++)
        wholeExtent_[i] = extent[i];
}

void ImageSource::setOrigin(const double * const origin)
{
    for (int i = 0; i < 3; i++)
        origin_[i] = origin[i];
}

void ImageSource::setSpacing(const double * const spacing)
{
    for (int i = 0; i < 3; i++)
        spacing_[i] = spacing[i];
}

int ImageSource::RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent_, 6);
    outInfo->Set(vtkDataObject::ORIGIN(), origin_, 3);
    outInfo->Set(vtkDataObject::SPACING(), spacing_, 3);

    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);

    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 1);

    return 1;
}

void ImageSource::ExecuteDataWithInformation(vtkDataObject *, vtkInformation *outInfo)
{
    int *execExt = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());

    vtkImageData *data = vtkImageData::GetData(outInfo);
    this->AllocateOutputData(data, outInfo, execExt);

    if (data->GetNumberOfPoints() <= 0)
      return;

    data->SetSpacing(spacing_[0], spacing_[1], spacing_[2]);

    int *outExt = data->GetExtent();

    data->GetPointData()->GetScalars()->SetName("myPointData");

    int bounds[3] = {outExt[1] - outExt[0],
                     outExt[3] - outExt[2],
                     outExt[5] - outExt[4]};

    vtkIdType outInc[3] = {1, 1, 1};
    data->GetContinuousIncrements(outExt, outInc[0], outInc[1], outInc[2]);

    float *outPtr = static_cast<float *>(data->GetScalarPointer(outExt[0], outExt[2], outExt[4]));

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    int offset = me ? bounds[2] : 0;
    for (int k = 0; k <= bounds[2]; k++)
    {
        for (int j = 0; j <= bounds[1]; j++)
        {
            for (int i = 0; i <= bounds[0]; i++)
            {
                *outPtr = (float)(i + j + k + offset);
                outPtr++;
            }
            outPtr += outInc[1];
        }
        outPtr += outInc[2];
    }

    vtkSmartPointer<vtkDoubleArray> ca = vtkSmartPointer<vtkDoubleArray>::New();
    ca->SetName("myCellData");
    double value = 1.0;
    for (int k = 0; k < bounds[2]; k++)
    {
        for (int j = 0; j < bounds[1]; j++)
        {
            for (int i = 0; i < bounds[0]; i++)
            {
                value = (double)(i + j + k + offset);
                ca->InsertNextTupleValue(&value);
            }
        }
    }
    data->GetCellData()->AddArray(ca);

}

int ImageSource::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
        return vtkImageAlgorithm::FillOutputPortInformation(port, info);

    return 0;
}

int main(int narg, char ** arg)
{
    MPI_Init(&narg, &arg);
    int me, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 2)
    {
        printf("Example runs only with either 1 or 2 processes\n");
        return -1;
    }

    vtkSmartPointer<vtkMPIController> controller = static_cast<vtkMPIController*>(vtkMultiProcessController::GetGlobalController());
    if (!controller)
    {
        controller = vtkMPIController::New();
        controller->Initialize();
        vtkMultiProcessController::SetGlobalController(controller);
    }

    vtkSmartPointer<ImageSource> source = new ImageSource;
    int zlow = (nprocs == 2 && me == 1) ? 1 : -10;
    int zhigh = (nprocs == 2 && me == 0) ? 0 : 10;
    int extent[6] = {-10, 10, -10, 10, zlow, zhigh};
    int whole_extent[6] = {-10, 10, -10, 10, -10, 10};
    double spacing[3] = {0.5, 1.0, 1.0};
    double origin[3] = {5., 0., 0.};
    source->setWholeExtent(whole_extent);
    source->setSpacing(spacing);
    source->setOrigin(origin);
    source->SetUpdateExtent(extent);
    vtkSmartPointer<vtkAlgorithmOutput> idport = source->GetOutputPort();

    vtkSmartPointer<vtkXMLPImageDataWriter> writer = vtkSmartPointer<vtkXMLPImageDataWriter>::New();
    writer->SetFileName("test.pvti");
    writer->SetInputConnection(idport);
    writer->SetNumberOfPieces(nprocs);
    writer->SetStartPiece(me);
    writer->SetEndPiece(me);
    writer->SetWriteSummaryFile(me == 0 ? 1 : 0);
    writer->Write();

    MPI_Finalize();

    return 0;
}
