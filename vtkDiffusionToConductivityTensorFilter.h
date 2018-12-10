// .NAME vtkDiffusionToConductivityTensorFilter - convert diff.-tensors to cond.-tensors

// .SECTION Description
// A Class to convert diffusion tensors to conductivity tensors.

// .SECTION Thanks
// noone

#ifndef vtkDiffusionToConductivityTensorFilter_h
#define vtkDiffusionToConductivityTensorFilter_h

#include "vtkDataSetAlgorithm.h"

class vtkDiffusionToConductivityTensorFilter : public vtkDataSetAlgorithm
{
public:
  vtkGetMacro(conductivity, double);
  vtkSetMacro(conductivity, double);
  vtkGetMacro(ratio, double);
  vtkSetMacro(ratio, double);

  vtkTypeMacro(vtkDiffusionToConductivityTensorFilter, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkDiffusionToConductivityTensorFilter *New();

  const double THRESHOLD = 1e-50;

protected:
  vtkDiffusionToConductivityTensorFilter();
  ~vtkDiffusionToConductivityTensorFilter();

  /* implementation of algorithm */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  void SetActiveTensors(int, int, int, int, const char *);

protected:
  double conductivity;
  double ratio;

private:
  vtkDiffusionToConductivityTensorFilter(const vtkDiffusionToConductivityTensorFilter&);  // Not implemented.
  void operator=(const vtkDiffusionToConductivityTensorFilter&);  // Not implemented.
};

#endif
