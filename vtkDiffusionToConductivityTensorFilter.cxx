#include "vtkDiffusionToConductivityTensorFilter.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkAppendPolyData.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"

#include "vtkCompositeDataIterator.h"
#include "vtkSelectEnclosedPoints.h"
#include "vtkSmartPointer.h"

#include <limits>
#include <csignal>

vtkStandardNewMacro(vtkDiffusionToConductivityTensorFilter);

vtkDiffusionToConductivityTensorFilter::vtkDiffusionToConductivityTensorFilter()
{
    this->SetNumberOfInputPorts(1);

    // by default, process active point tensors
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);
}

vtkDiffusionToConductivityTensorFilter::~vtkDiffusionToConductivityTensorFilter()
{
}

int vtkDiffusionToConductivityTensorFilter::FillInputPortInformation
(
  int               port,
  vtkInformation*   info
)
{
    info->Remove( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() );
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );

    return 1;
}


int vtkDiffusionToConductivityTensorFilter::RequestData
(
  vtkInformation        *vtkNotUsed(request),
  vtkInformationVector  **inputVector,
  vtkInformationVector  *outputVector
  )
{
    // get the info objects
    vtkInformation *inInfoP0 = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo  = outputVector->GetInformationObject(0);

    // get the input and output
    vtkDataSet *tensors = vtkDataSet::SafeDownCast( inInfoP0->Get( vtkDataObject::DATA_OBJECT() ) );
    vtkDataSet *output  = vtkDataSet::SafeDownCast( outInfo->Get( vtkDataObject::DATA_OBJECT() ) );

    vtkPointData *tensorsPD = tensors->GetPointData();
    vtkPointData *outputPD  = output->GetPointData();

    if ( !tensorsPD )
    {
        vtkErrorMacro(<<"No input tensor-data!");
        return 0;
    }

    // First, copy the structure of the input to the output
    // we will assing the field-values later
    output->CopyStructure( tensors );

    // proceed only, if both mask and tensor-field have data-points (attributes)
    if( tensorsPD->GetNumberOfArrays() > 0)
    {
        // retrieve the data-points
        vtkDataArray *tensorData  = tensorsPD->GetArray(0);
        vtkDataArray *outputData  = vtkDataArray::CreateDataArray(VTK_DOUBLE);

        vtkIdType numOfTuples  = tensorsPD->GetNumberOfTuples();
        int componentsPerValue = tensorsPD->GetNumberOfComponents();
        outputData->SetNumberOfComponents( 9 );         // from doc: "Make sure that this is set before allocation. "
        outputData->SetNumberOfTuples( numOfTuples );   // allocation
        outputData->SetName( "Conductivity-Tensors" );

        // some tensor-indepentant calculations, we can solve before the
        // actual processing:
        // - calculate the eigenvalues of the conductivity-tensor
        //   according to wang's constraint
        //      isotropicConductivity^2 = longitudinal * transversal
        //           with
        //   transversal : longitudinal =  RATIO : 1
        // - create diagonal matrix using these eigenvalues
        double isoCondSqr   = conductivity * conductivity;
        double transversal  = sqrt( isoCondSqr * ratio );    // auxiliary-direction
        double longitudinal = isoCondSqr / transversal;      // main-direction

        double kroneckerScaled[3][3]
            = {{longitudinal,0, 0}, {0,transversal,0}, {0,0,transversal}};

        // prepare variables; we have to use different formats for
        // the various math-functions we call... -.-
            //tensor
        double tensor[9];                   // this is how we get the data from the field
        double zeroes[9] = {0,0,0,0,0,0,0,0,0};
        double (*tensor9x1)[9]= nullptr;    // this is for writing the data back to the field (convert 3x3 to 9x1 matrix)
        double *tensorMat[3];               // this is the format vtKMath:Jacobi requires
        double m0[3];double m1[3];double m2[3];
        tensorMat[0] = m0;
        tensorMat[1] = m1;
        tensorMat[2] = m2;
            // eigenvalue, eigenvector
        double eigenVectors3x3[3][3];   // this is the format vtkMath::Multiply3x3 requires
        double eigenVectors3x3_T[3][3]; // -"-
        double temp[3][3];              // temporary variable for results of vtkMath::Multiply3x3
        double condTensor[3][3];        // the result, the conductivity tensor
        double eigenValues[3];          // this is the format vtkMath::Jacobi requires
        double *eigenVectors[3];        // -"-
        eigenVectors[0] = eigenVectors3x3[0];   // assign every row to the pointer-array
        eigenVectors[1] = eigenVectors3x3[1];
        eigenVectors[2] = eigenVectors3x3[2];

        unsigned char zeroCtr = 0;

        for( vtkIdType t = 0; t < numOfTuples; t++)
        {
            // retrieve diffusion tensor
            tensorData->GetTuple( t, tensor );

                // copy tensor since it will be modified by Jacobi
                // and is expected as 2D array
            for( unsigned char i = 0; i < 3; i++ )
            {
                for( unsigned char j = 0; j < 3; j++)
                {
                    if( fabs(tensor[3*i+j]) < THRESHOLD )
                    {
                        zeroCtr++;
                        tensorMat[i][j] = 0;
                    }
                    else
                    {
                        tensorMat[i][j] = tensor[3*i+j];
                    }
                }
            }

            if( zeroCtr < 9 )
            {
                // eigenvector & -value decomposition
                // from documentation:
                // "Resulting eigenvalues/vectors are sorted in decreasing
                // order; eigenvectors are normalized."
                // => exactly what we need!
                vtkMath::Jacobi( tensorMat, eigenValues, eigenVectors );

                // we transform the cond-tensor eigenvalue matrix using
                // the eigenvectors of the diff-tensor:
                //  cond_tensor = eVec_D * eValues_C * eVec_D^T
                vtkMath::Transpose3x3( eigenVectors3x3, eigenVectors3x3_T );
                vtkMath::Multiply3x3( eigenVectors3x3, kroneckerScaled, temp );
                vtkMath::Multiply3x3( temp, eigenVectors3x3_T, condTensor );

                tensor9x1 = reinterpret_cast<double(*)[9]>(condTensor);

                outputData->SetTuple( t, *tensor9x1 );
            }
            else
            {
                outputData->SetTuple( t, zeroes );
            }

            zeroCtr = 0;
        }

        outputPD->AddArray( outputData );
    }
    else
    {
        vtkErrorMacro(<<"The input fields does not conatin any data!");
    }

  return 1;
}

void vtkDiffusionToConductivityTensorFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
