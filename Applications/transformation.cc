/* The Image Registration Toolkit (IRTK)
 *
 * Copyright 2008-2015 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. */

#include <irtkImage.h>
#include <irtkFileToImage.h>
#include <irtkRegionFilter.h>
#include <irtkTransformation.h>
#include <memory>


// ===========================================================================
// Help
// ===========================================================================

void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <source> <output> [options]" << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -dofin <file>    Transformation or 'Id'/'Identity'" << endl;
  cout << "  -target <file>   Target image" << endl;
  cout << "  -2d              Project transformation in 2D" << endl;
  cout << "  -Rx1 <value>     Region of interest" << endl;
  cout << "  -Ry1 <value>     Region of interest" << endl;
  cout << "  -Rz1 <value>     Region of interest" << endl;
  cout << "  -Rx2 <value>     Region of interest" << endl;
  cout << "  -Ry2 <value>     Region of interest" << endl;
  cout << "  -Rz2 <value>     Region of interest" << endl;
  cout << "  -Tx1 <value>     Region of interest in target image" << endl;
  cout << "  -Ty1 <value>     Region of interest in target image" << endl;
  cout << "  -Tz1 <value>     Region of interest in target image" << endl;
  cout << "  -Tx2 <value>     Region of interest in target image" << endl;
  cout << "  -Ty2 <value>     Region of interest in target image" << endl;
  cout << "  -Tz2 <value>     Region of interest in target image" << endl;
  cout << "  -Sx1 <value>     Region of interest in source image" << endl;
  cout << "  -Sy1 <value>     Region of interest in source image" << endl;
  cout << "  -Sz1 <value>     Region of interest in source image" << endl;
  cout << "  -Sx2 <value>     Region of interest in source image" << endl;
  cout << "  -Sy2 <value>     Region of interest in source image" << endl;
  cout << "  -Sz2 <value>     Region of interest in source image" << endl;
  cout << "  -Tp <value>      Target padding value" << endl;
  cout << "  -Sp <value>      Source padding value" << endl;
  cout << "  -invert          Invert transformation" << endl;
  cout << "  -nn              Nearst Neighbor interpolation" << endl;
  cout << "  -linear          Linear interpolation" << endl;
  cout << "  -bspline         B-spline interpolation" << endl;
  cout << "  -cspline         Cubic spline interpolation" << endl;
  cout << "  -sinc            Sinc interpolation" << endl;
  cout << "  -sbased          Shape based interpolation" << endl;
  cout << "  -matchInputType  Make the output data type (short, float, etc.)" << endl;
  cout << "                   the same as that of the input image regardless" << endl;
  cout << "                   of the data type of the target (if specified)." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// ===========================================================================
// Main
// ===========================================================================

int main(int argc, char **argv)
{
  // Parse positional arguments
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  // Read image
  irtkFileToImage *inputReader = irtkFileToImage::New(input_name);
  int sourceType = inputReader->GetDataType();
  std::unique_ptr<irtkImage> source(inputReader->GetOutput());
  delete inputReader;

  // Parse optional arguments
  const char           *dof_name      = NULL;
  bool                  dof_invert    = false;
  const char           *dofin_name    = NULL;
  const char           *target_name   = NULL;
  irtkInterpolationMode interpolation = Interpolation_NN;

  int target_x1 = -1;
  int target_y1 = -1;
  int target_z1 = -1;
  int target_x2 = -1;
  int target_y2 = -1;
  int target_z2 = -1;
  int source_x1 = -1;
  int source_y1 = -1;
  int source_z1 = -1;
  int source_x2 = -1;
  int source_y2 = -1;
  int source_z2 = -1;

  double target_t = numeric_limits<double>::quiet_NaN();
  double source_t = numeric_limits<double>::quiet_NaN();

  int  source_padding  = 0;
  int  target_padding  = MIN_GREY;
  bool matchSourceType = false;
  bool invert          = false;
  bool twod            = false;

  for (ALL_OPTIONS) {
      if      (OPTION("-dof")    ) dof_name       = ARGUMENT, dof_invert = false;
      else if (OPTION("-dof_i")  ) dof_name       = ARGUMENT, dof_invert = true;
      else if (OPTION("-dofin")  ) dofin_name     = ARGUMENT;
      else if (OPTION("-target") ) target_name    = ARGUMENT;
      else if (OPTION("-Rx1")    ) target_x1 = source_x1 = atoi(ARGUMENT);
      else if (OPTION("-Rx2")    ) target_x2 = source_x2 = atoi(ARGUMENT);
      else if (OPTION("-Ry1")    ) target_y1 = source_y1 = atoi(ARGUMENT);
      else if (OPTION("-Ry2")    ) target_y2 = source_y2 = atoi(ARGUMENT);
      else if (OPTION("-Rz1")    ) target_z1 = source_z1 = atoi(ARGUMENT);
      else if (OPTION("-Rz2")    ) target_z2 = source_z2 = atoi(ARGUMENT);
      else if (OPTION("-Tx1")    ) target_x1      = atoi(ARGUMENT);
      else if (OPTION("-Tx2")    ) target_x2      = atoi(ARGUMENT);
      else if (OPTION("-Ty1")    ) target_y1      = atoi(ARGUMENT);
      else if (OPTION("-Ty2")    ) target_y2      = atoi(ARGUMENT);
      else if (OPTION("-Tz1")    ) target_z1      = atoi(ARGUMENT);
      else if (OPTION("-Tz2")    ) target_z2      = atoi(ARGUMENT);
      else if (OPTION("-Sx1")    ) source_x1      = atoi(ARGUMENT);
      else if (OPTION("-Sx2")    ) source_x2      = atoi(ARGUMENT);
      else if (OPTION("-Sy1")    ) source_y1      = atoi(ARGUMENT);
      else if (OPTION("-Sy2")    ) source_y2      = atoi(ARGUMENT);
      else if (OPTION("-Sz1")    ) source_z1      = atoi(ARGUMENT);
      else if (OPTION("-Sz2")    ) source_z2      = atoi(ARGUMENT);
      else if (OPTION("-Tt")     ) target_t       = atof(ARGUMENT);
      else if (OPTION("-St")     ) source_t       = atof(ARGUMENT);
      else if (OPTION("-Tp")     ) target_padding = atoi(ARGUMENT);
      else if (OPTION("-Sp")     ) source_padding = atoi(ARGUMENT);
      else if (OPTION("-invert") ) invert         = true;
      else if (OPTION("-2d")     ) twod           = true;
      else if (OPTION("-nn")     ) interpolation  = Interpolation_NN;
      else if (OPTION("-linear") ) interpolation  = Interpolation_Linear;
      else if (OPTION("-bspline")) interpolation  = Interpolation_BSpline;
      else if (OPTION("-cspline")) interpolation  = Interpolation_CSpline;
      else if (OPTION("-sinc")   ) interpolation  = Interpolation_Sinc;
      else if (OPTION("-sbased") ) interpolation  = Interpolation_SBased;
      else if (OPTION("-matchInputType")) matchSourceType = true;
      else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Apply affine header transformation
  if (dof_name) {
    std::unique_ptr<irtkTransformation> t(irtkTransformation::New(dof_name));
    irtkHomogeneousTransformation *lin = dynamic_cast<irtkHomogeneousTransformation *>(t.get());
    if (lin) {
      irtkMatrix mat = lin->GetMatrix();
      if (dof_invert) mat.Invert();
      source->PutAffineMatrix(mat, true);
    } else {
      cerr << "Error: -dof transformation can only be affine" << endl;
      exit(1);
    }
  }

  // Read target image
  std::unique_ptr<irtkImage> target;
  int targetType = -1;
  if (target_name) {
    irtkFileToImage *targetReader = irtkFileToImage::New(target_name);
    targetType = targetReader->GetDataType();
    target.reset(targetReader->GetOutput());
    delete targetReader;
  }

  // Instantiate interpolator
  std::unique_ptr<irtkInterpolateImageFunction> interpolator(irtkInterpolateImageFunction::New(interpolation));

  // Initialize output image
  if (target.get()) {
    // Change type of target image or warn about differing output type
    if (sourceType != targetType) {
      if (matchSourceType) {
        std::unique_ptr<irtkImage> tmp(irtkImage::New(sourceType));
        *tmp = *target;
        target.reset(tmp.release());
        targetType = sourceType;
      } else {
        cerr << "\nWarning! source and target image have different data types:" << endl;
        cerr << "Converting data from " << DataTypeName(sourceType) << " to " << DataTypeName(targetType) << ".\n" << endl;
      }
    }
    // Ensure that output image has same number of channels/frames than input
    if (target->T() != source->T()) {
      std::unique_ptr<irtkImage> tmp(irtkImage::New(targetType));
      tmp->Initialize(target->Attributes(), source->T());
      tmp->PutTSize(source->GetTSize());
      for (int l = 0; l < tmp->T(); ++l)
      for (int k = 0; k < tmp->Z(); ++k)
      for (int j = 0; j < tmp->Y(); ++j)
      for (int i = 0; i < tmp->X(); ++i) {
        tmp->PutAsDouble(i, j, k, l, target->GetAsDouble(i, j, k, 0));
      }
      target.reset(tmp.release());
    }
  // If there is no target image just copy the source image
  } else {
    target.reset(irtkImage::New(source.get()));
  }

  // Set temporal offset
  if (!IsNaN(target_t)) target->PutTOrigin(target_t);
  if (!IsNaN(source_t)) source->PutTOrigin(source_t);

  // Compute region of interest for target image
  if (target_x1 == -1) target_x1 = 0;
  if (target_y1 == -1) target_y1 = 0;
  if (target_z1 == -1) target_z1 = 0;
  if (target_x2 == -1) target_x2 = target->GetX();
  if (target_y2 == -1) target_y2 = target->GetY();
  if (target_z2 == -1) target_z2 = target->GetZ();

  // Compute region of interest for source image
  if (source_x1 == -1) source_x1 = 0;
  if (source_y1 == -1) source_y1 = 0;
  if (source_z1 == -1) source_z1 = 0;
  if (source_x2 == -1) source_x2 = source->GetX();
  if (source_y2 == -1) source_y2 = source->GetY();
  if (source_z2 == -1) source_z2 = source->GetZ();

  // If there is an region of interest for the target image, use it
  if ((target_x1 != 0) || (target_x2 != target->GetX()) ||
      (target_y1 != 0) || (target_y2 != target->GetY()) ||
      (target_z1 != 0) || (target_z2 != target->GetZ())) {
    irtkRegionFilter *region = new irtkRegionFilter;
    region->SetInput (target.get());
    region->SetOutput(target.get());
    region->PutRegion(target_x1, target_y1, target_z1, 0, target_x2, target_y2, target_z2, target->GetT());
    region->Run();
    delete region;
  }

  // If there is an region of interest for the source image, use it
  if ((source_x1 != 0) || (source_x2 != source->GetX()) ||
      (source_y1 != 0) || (source_y2 != source->GetY()) ||
      (source_z1 != 0) || (source_z2 != source->GetZ())) {
    irtkRegionFilter *region = new irtkRegionFilter;
    region->SetInput (source.get());
    region->SetOutput(source.get());
    region->PutRegion(source_x1, source_y1, source_z1, 0, source_x2, source_y2, source_z2, source->GetT());
    region->Run();
    delete region;
  }

  // Instantiate image transformation
  std::unique_ptr<irtkTransformation> transformation;
  if (dofin_name == NULL || strcmp(dofin_name, "identity") == 0
                         || strcmp(dofin_name, "Identity") == 0
                         || strcmp(dofin_name, "Id")       == 0) {
    // Create identity transformation
    transformation.reset(new irtkRigidTransformation);
  } else {
    // Read transformation
    transformation.reset(irtkTransformation::New(dofin_name));
  }

  // Create image transformation filter
  std::unique_ptr<irtkImageTransformation> imagetransformation(irtkImageTransformation::New(transformation.get()));

  imagetransformation->SetInput(source.get(), transformation.get());
  imagetransformation->SetOutput(target.get());
  imagetransformation->PutTargetPaddingValue(target_padding);
  imagetransformation->PutSourcePaddingValue(source_padding);
  imagetransformation->PutInterpolator(interpolator.get());

  if (invert) imagetransformation->InvertOn();
  else        imagetransformation->InvertOff();

  if (twod)   imagetransformation->TwoDOn();
  else        imagetransformation->TwoDOff();

  // Transform source image
  imagetransformation->Run();
  if (imagetransformation->NumberOfSingularPoints() > 0) {
    cerr << "Warning: Transformation is non-invertible at "
         << imagetransformation->NumberOfSingularPoints() << " point";
    if (imagetransformation->NumberOfSingularPoints() > 1) cerr << 's';
    cerr << endl;
  }

  // Write the final transformation estimate
  target->Write(output_name);

  return 0;
}
