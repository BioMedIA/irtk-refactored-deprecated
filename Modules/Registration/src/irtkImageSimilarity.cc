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

#include <irtkImageSimilarity.h>
#include <irtkVoxelFunction.h>


// =============================================================================
// Factory
// =============================================================================

#include <irtkSumOfSquaredIntensityDifferences.h>
#include <irtkIntensityCrossCorrelation.h>
#include <irtkNormalizedIntensityCrossCorrelation.h>
//#include <irtkJointImageEntropy.h>
#include <irtkMutualImageInformation.h>
#include <irtkNormalizedMutualImageInformation.h>
#include <irtkCosineOfNormalizedGradientField.h>

// -----------------------------------------------------------------------------
irtkImageSimilarity *irtkImageSimilarity::New(irtkSimilarityMeasure metric)
{
  switch (metric) {
    case SSD:      return new irtkSumOfSquaredIntensityDifferences();
    case  CC:      return new irtkIntensityCrossCorrelation();
    case NCC:      return new irtkNormalizedIntensityCrossCorrelation();
//    case  JE: return new irtkJointImageEntropy();
    case  MI:      return new irtkMutualImageInformation();
    case NMI:      return new irtkNormalizedMutualImageInformation();
    case NGF_COS:  return new irtkCosineOfNormalizedGradientField();
    default:
      cerr << "irtkImageSimilarity::New: Unknown image (dis-)similarity measure: " << metric << endl;
      exit(1);
  }
  return NULL;
}

// =============================================================================
// Auxiliary functor classes for parallel execution
// =============================================================================

namespace irtkImageSimilarityUtils {


// -----------------------------------------------------------------------------
/// Voxel function used by MultiplyByImageGradient to post-multiply similarity
/// gradient by transformed image gradient according to the chain rule.
class MultiplySimilarityGradientByImageGradient : public irtkVoxelFunction
{
private:

  int _NumberOfVoxels;
  int _dx, _dy, _dz;

public:

  MultiplySimilarityGradientByImageGradient(const irtkRegisteredImage *image)
  :
    _NumberOfVoxels(image->NumberOfVoxels()),
    _dx            (image->Offset(irtkRegisteredImage::Dx)),
    _dy            (image->Offset(irtkRegisteredImage::Dy)),
    _dz            (image->Offset(irtkRegisteredImage::Dz))
  {}

  template <class Image, class TScalar, class TReal>
  void operator ()(const Image &, int, const TScalar *dI, TReal *gradient)
  {
    (*gradient) *= dI[_dx]; gradient += _NumberOfVoxels;
    (*gradient) *= dI[_dy]; gradient += _NumberOfVoxels;
    (*gradient) *= dI[_dz];
  }
};

// -----------------------------------------------------------------------------
/// Determine maximum norm of voxel-wise image similarity gradient
class MaxVoxelWiseSimilarityGradient : public irtkVoxelReduction
{
private:

  double    _norm;
  const int _x, _y, _z;

public:

  /// Constructor
  MaxVoxelWiseSimilarityGradient(const irtkBaseImage *img)
  :
    _norm(.0), _x(0), _y(img->X() * img->Y() * img->Z()), _z(2 * _y)
  {}

  /// Split "constructor"
  void split(const MaxVoxelWiseSimilarityGradient &other)
  {
    _norm = other._norm;
  }

  /// Join results
  void join(const MaxVoxelWiseSimilarityGradient &other)
  {
    if (other._norm > _norm) _norm = other._norm;
  }

  /// Get maximum norm
  double Norm() const { return sqrt(_norm); }

  /// Compute norm of similarity gradient
  template <class Image, class TReal>
  void operator ()(const Image &, int, TReal *gradient)
  {
    double norm = pow(gradient[_x], 2) + pow(gradient[_y], 2) + pow(gradient[_z], 2);
    if (norm > _norm) _norm = norm;
  }
};

// -----------------------------------------------------------------------------
/// Determine maximum norm of node-based image similarity gradient
class MaxNodeBasedSimilarityGradient
{
private:

  const irtkFreeFormTransformation *_FFD;
  const double                     *_Gradient;
  double                            _MaxNorm;

public:

  /// Constructor
  MaxNodeBasedSimilarityGradient(const irtkFreeFormTransformation *ffd,
                                 const double                     *gradient)
  :
    _FFD(ffd), _Gradient(gradient), _MaxNorm(.0)
  {}

  /// Copy constructor
  MaxNodeBasedSimilarityGradient(const MaxNodeBasedSimilarityGradient &other)
  :
    _FFD     (other._FFD),
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm)
  {}

  /// Split constructor
  MaxNodeBasedSimilarityGradient(const MaxNodeBasedSimilarityGradient &other, split)
  :
    _FFD     (other._FFD),
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm)
  {}

  /// Join results
  void join(const MaxNodeBasedSimilarityGradient &other)
  {
    if (other._MaxNorm > _MaxNorm) _MaxNorm = other._MaxNorm;
  }

  /// Maximum norm
  double Norm() const { return sqrt(_MaxNorm); }

  /// Determine maximum norm of specified control point gradients
  void operator()(const blocked_range<int> &re)
  {
    double norm;
    int    x, y, z;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToDOFs(cp, x, y, z);
      norm = pow(_Gradient[x], 2) + pow(_Gradient[y], 2) + pow(_Gradient[z], 2);
      if (norm > _MaxNorm) _MaxNorm = norm;
    }
  }
};

// -----------------------------------------------------------------------------
/// Normalize voxel-wise image similarity gradient
class NormalizeVoxelWiseSimilarityGradient : public irtkVoxelFunction
{
private:

  double    _Sigma;
  const int _x, _y, _z;

public:

  /// Constructor
  NormalizeVoxelWiseSimilarityGradient(const irtkBaseImage *img, double sigma = .0)
  :
    _Sigma(sigma), _x(0), _y(img->X() * img->Y() * img->Z()), _z(2 * _y)
  {}

  /// Normalize similarity gradient
  template <class Image, class TReal>
  void operator ()(const Image &, int, TReal *gradient)
  {
    double norm = pow(gradient[_x], 2) + pow(gradient[_y], 2) + pow(gradient[_z], 2);
    if (norm) {
      norm = sqrt(norm) + _Sigma;
      gradient[_x] /= norm;
      gradient[_y] /= norm;
      gradient[_z] /= norm;
    }
  }
};

// -----------------------------------------------------------------------------
/// Normalize node-based image similarity gradient
class NormalizeNodeBasedSimilarityGradient
{
private:

  const irtkFreeFormTransformation *_FFD;
  double                           *_Gradient;
  double                            _Sigma;

public:

  /// Constructor
  NormalizeNodeBasedSimilarityGradient(const irtkFreeFormTransformation *ffd,
                                       double                           *gradient,
                                       double                            sigma)
  :
    _FFD(ffd), _Gradient(gradient), _Sigma(sigma)
  {}

  /// Copy constructor
  NormalizeNodeBasedSimilarityGradient(const NormalizeNodeBasedSimilarityGradient &other)
  :
    _FFD     (other._FFD),
    _Gradient(other._Gradient),
    _Sigma   (other._Sigma)
  {}

  /// Normalize node-based similarity gradient
  void operator ()(const blocked_range<int> &re) const
  {
    double norm;
    int    x, y, z;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToDOFs(cp, x, y, z);
      norm = pow(_Gradient[x], 2) + pow(_Gradient[y], 2) + pow(_Gradient[z], 2);
      if (norm) {
        norm = sqrt(norm) + _Sigma;
        _Gradient[x] /= norm;
        _Gradient[y] /= norm;
        _Gradient[z] /= norm;
      }
    }
  }
};


} // namespace irtkImageSimilarityUtils
using namespace irtkImageSimilarityUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkImageSimilarity::irtkImageSimilarity(const char *name, double weight)
:
  irtkDataFidelity(name, weight),
  _Target                  (new irtkRegisteredImage()),
  _Source                  (new irtkRegisteredImage()),
  _GradientWrtTarget       (NULL),
  _GradientWrtSource       (NULL),
  _Gradient                (NULL),
  _NumberOfVoxels          (0),
  _UseApproximateGradient  (false),
  _VoxelWisePreconditioning(.0),
  _NodeBasedPreconditioning(.0),
  _InitialUpdate           (false)
{
  _ParameterPrefix.push_back("Image (dis-)similarity ");
  _ParameterPrefix.push_back("Image dissimilarity ");
  _ParameterPrefix.push_back("Image similarity ");
  _ParameterPrefix.push_back("(Dis-)similarity ");
  _ParameterPrefix.push_back("Dissimilarity ");
  _ParameterPrefix.push_back("Similarity ");
}

// -----------------------------------------------------------------------------
irtkImageSimilarity::irtkImageSimilarity(const irtkImageSimilarity &other)
:
  irtkDataFidelity(other),
  _Target           (other._Target ? new irtkRegisteredImage(*other._Target) : NULL),
  _Source           (other._Source ? new irtkRegisteredImage(*other._Source) : NULL),
  _GradientWrtTarget       (NULL),
  _GradientWrtSource       (NULL),
  _Gradient                (NULL),
  _NumberOfVoxels          (other._NumberOfVoxels),
  _UseApproximateGradient  (other._UseApproximateGradient),
  _VoxelWisePreconditioning(other._VoxelWisePreconditioning),
  _NodeBasedPreconditioning(other._NodeBasedPreconditioning),
  _InitialUpdate           (other._InitialUpdate)
{
}

// -----------------------------------------------------------------------------
irtkImageSimilarity &irtkImageSimilarity::operator =(const irtkImageSimilarity &other)
{
  irtkDataFidelity::operator =(other);
  Delete(_Target);
  Delete(_Source);
  Delete(_GradientWrtTarget);
  Delete(_GradientWrtSource);
  Deallocate(_Gradient);
  _Target = other._Target ? new irtkRegisteredImage(*other._Target) : NULL;
  _Source = other._Source ? new irtkRegisteredImage(*other._Source) : NULL;
  _NumberOfVoxels           = other._NumberOfVoxels;
  _UseApproximateGradient   = other._UseApproximateGradient;
  _VoxelWisePreconditioning = other._VoxelWisePreconditioning;
  _NodeBasedPreconditioning = other._NodeBasedPreconditioning;
  _InitialUpdate            = other._InitialUpdate;
  return *this;
}

// -----------------------------------------------------------------------------
irtkImageSimilarity::~irtkImageSimilarity()
{
  Delete(_Target);
  Delete(_Source);
  Delete(_GradientWrtTarget);
  Delete(_GradientWrtSource);
  Deallocate(_Gradient);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageSimilarity::InitializeInput(const irtkImageAttributes &domain)
{
  _Target->Initialize(domain, _Target->Transformation() ? 4 : 1);
  _Source->Initialize(domain, _Source->Transformation() ? 4 : 1);
  _Domain = domain;
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::Initialize()
{
  // Free previously allocated memory
  Delete(_GradientWrtTarget);
  Delete(_GradientWrtSource);
  Deallocate(_Gradient);
  // Check if all inputs are set
  if (!_Target) {
    cerr << "irtkImageSimilarity::Initialize: Missing target image" << endl;
    exit(1);
  }
  if (!_Source) {
    cerr << "irtkImageSimilarity::Initialize: Missing source image" << endl;
    exit(1);
  }
  if (( _Mask &&  _Mask->IsEmpty()) ||
      (!_Mask && (_Domain._x == 0 || _Domain._y == 0 || _Domain._z == 0))) {
    cerr << "irtkImageSimilarity::Initialize: No image domain specified" << endl;
    exit(1);
  }
  // Initialize base class
  irtkDataFidelity::Initialize();
  // Initialize registered images
  this->InitializeInput(_Mask ? _Mask->Attributes() : _Domain);
  _InitialUpdate = true; // i.e., initialize image content upon first Update
  // Allocate memory for temporary similarity gradient
  if (_NodeBasedPreconditioning > .0) {
    const irtkTransformation *T1 = _Target->Transformation();
    const irtkTransformation *T2 = _Source->Transformation();
    const int n = max(T1 ? T1->NumberOfDOFs() : 0,
                      T2 ? T2->NumberOfDOFs() : 0);
    Allocate(_Gradient, n);
  }
  // Number of voxels per registered image
  _NumberOfVoxels = _Domain.NumberOfSpatialPoints();
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkImageSimilarity::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);

  if (name == "Approximate gradient") {
    return FromString(value, _UseApproximateGradient);
  }
  if (name == "Preconditioning (voxel-wise)") {
    return FromString(value, _VoxelWisePreconditioning);
  }
  if (name == "Preconditioning (node-based)" ||
      name == "Preconditioning") {
    return FromString(value, _NodeBasedPreconditioning);
  }
  if (name == "Blurring of 1st order image derivatives" ||
      name == "Blurring of image jacobian" ||
      name == "Blurring of image gradient") {
    double sigma;
    if (!FromString(value, sigma)) return false;
    _Target->GradientSigma(sigma);
    _Source->GradientSigma(sigma);
    return true;
  }
  if (name == "Blurring of 2nd order image derivatives" ||
      name == "Blurring of image hessian") {
    double sigma;
    if (!FromString(value, sigma)) return false;
    _Target->HessianSigma(sigma);
    _Source->HessianSigma(sigma);
    return true;
  }

  return irtkDataFidelity::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkImageSimilarity::Parameter() const
{
  irtkParameterList params = irtkDataFidelity::Parameter();
  if (!_Name.empty()) {
    Insert(params, _Name + " approximate gradient",         ToString(_UseApproximateGradient));
    Insert(params, _Name + " preconditioning (voxel-wise)", ToString(_VoxelWisePreconditioning));
    Insert(params, _Name + " preconditioning (node-based)", ToString(_NodeBasedPreconditioning));
    Insert(params, _Name + " blurring of image gradient",   ToString(_Target->GradientSigma()));
    Insert(params, _Name + " blurring of image hessian",    ToString(_Target->HessianSigma()));
  }
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageSimilarity::Update(bool gradient)
{
  if (_InitialUpdate || _Target->Transformation()) {
    _Target->Update(true, gradient, false, _InitialUpdate);
  }
  if (_InitialUpdate || _Source->Transformation()) {
    _Source->Update(true, gradient, false, _InitialUpdate);
  }
  _InitialUpdate = false;
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::Exclude(const blocked_range3d<int> &)
{
  // By default, call Update upon Include
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::Include(const blocked_range3d<int> &)
{
  this->Update(false);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::MultiplyByImageGradient(const irtkRegisteredImage *image,
                                                  GradientImageType         *gradient)
{
  // Copy (dSimilarity / dI) from x component also to y and z components
  const int nbytes = image->NumberOfVoxels() * sizeof(GradientType);
  GradientType *gx = gradient->GetPointerToVoxels();
  memcpy(gradient->GetPointerToVoxels(0, 0, 0, 1), gx, nbytes);
  memcpy(gradient->GetPointerToVoxels(0, 0, 0, 2), gx, nbytes);
  // Apply chain rule (dSimilarity / dy) = (dSimilarity / dI) * (dI / dy)
  // where y = T(x) to obtain the non-parametric similarity gradient
  MultiplySimilarityGradientByImageGradient times_dIdx(image);
  ParallelForEachVoxel(image, gradient, times_dIdx);
}

// -----------------------------------------------------------------------------
bool irtkImageSimilarity::NonParametricGradient(const irtkRegisteredImage *,
                                                GradientImageType         *)
{
  return false; // By default, ApproximateGradient using finite differences
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::NormalizeGradient(GradientImageType *gradient)
{
  IRTK_START_TIMING();

  // Determine maximum norm of control point gradients
  MaxVoxelWiseSimilarityGradient maximum(gradient);
  ParallelForEachVoxel(gradient, maximum);

  // Sigma value used to suppress noise
  const double sigma = _VoxelWisePreconditioning * maximum.Norm();

  // Normalize control point gradients to be possibly similar
  NormalizeVoxelWiseSimilarityGradient norm(gradient, sigma);
  ParallelForEachVoxel(gradient, norm);

  IRTK_DEBUG_TIMING(2, "normalization of voxel-wise (dis-)similarity gradient");
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::ParametricGradient(const irtkRegisteredImage *image,
                                             GradientImageType         *np_gradient,
                                             double                    *gradient,
                                             double                     weight)
{
  const irtkTransformation *T = image->Transformation();
  irtkAssert(T != NULL, "image is being transformed");
  const irtkWorldCoordsImage *i2w = image->ImageToWorld();
  const double                t0  = image->GetTOrigin();
  np_gradient->PutTOrigin(image->InputImage()->GetTOrigin());
  T->ParametricGradient(np_gradient, gradient, i2w, t0, weight);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::ApproximateGradient(irtkRegisteredImage        *image,
                                              irtkFreeFormTransformation *ffd,
                                              double *gradient, double step,
                                              double weight)
{
  weight /= 2.0 * step;
  double a, b, value;
  int i1, j1, k1, i2, j2, k2, dof[3];
  for (int cp = 0; cp < ffd->NumberOfCPs(); ++cp) {
    if (ffd->IsActive(cp) &&
        ffd->BoundingBox(image, cp, i1, j1, k1, i2, j2, k2)) {
      blocked_range3d<int> region(k1, k2+1, j1, j2+1, i1, i2+1);
      ffd->IndexToDOFs(cp, dof[0], dof[1], dof[2]);
      for (int i = 0; i < 3; ++i) {
        if (ffd->GetStatus(dof[i]) == Active) {
          value = ffd->Get(dof[i]);

          ffd->Put(dof[i], value + step);
          this->Exclude(region);
          image->Update(region);
          this->Include(region);
          a = this->Evaluate();

          ffd->Put(dof[i], value - step);
          this->Exclude(region);
          image->Update(region);
          this->Include(region);
          b = this->Evaluate();

          ffd->Put(dof[i], value);
          this->Exclude(region);
          image->Update(region);
          this->Include(region);

          gradient[dof[i]] += weight * (a - b);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::ApproximateGradient(irtkRegisteredImage *image,
                                              double *gradient, double step,
                                              double weight)
{
  IRTK_START_TIMING();

  irtkTransformation *T = const_cast<irtkTransformation *>(image->Transformation());
  irtkAssert(T != NULL, "image is being transformed");
  irtkFreeFormTransformation   *ffd  = NULL;
  irtkMultiLevelTransformation *mffd = NULL;
  (mffd = dynamic_cast<irtkMultiLevelTransformation *>(T)) ||
  (ffd  = dynamic_cast<irtkFreeFormTransformation   *>(T));

  if (mffd) {
    for (int i = 0; i < mffd->NumberOfLevels(); ++i) {
      if (mffd->LocalTransformationIsActive(i)) {
        ffd = mffd->GetLocalTransformation(i);
        this->ApproximateGradient(image, ffd, gradient, step, weight);
        gradient += ffd->NumberOfDOFs();
      }
    }
  } else if (ffd) {
    this->ApproximateGradient(image, ffd, gradient, step, weight);
  } else {
    weight /= 2.0 * step;
    double a, b, value;
    blocked_range3d<int> region(0, image->Z(), 0, image->Y(), 0, image->X());
    for (int dof = 0; dof < T->NumberOfDOFs(); ++dof) {
      if (T->GetStatus(dof) == Active) {
        value = T->Get(dof);

        T->Put(dof, value + step);
        this->Exclude(region);
        image->Update(region);
        this->Include(region);
        a = this->Evaluate();

        T->Put(dof, value - step);
        this->Exclude(region);
        image->Update(region);
        this->Include(region);
        b = this->Evaluate();

        T->Put(dof, value);
        this->Exclude(region);
        image->Update(region);
        this->Include(region);

        gradient[dof] += weight * (a - b);
      }
    }
  }

  IRTK_DEBUG_TIMING(2, "approximation of similarity gradient");
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::NormalizeGradient(const irtkRegisteredImage *image, double *gradient)
{
  const irtkMultiLevelTransformation *mffd = NULL;
  const irtkFreeFormTransformation   *affd = NULL;

  (mffd = dynamic_cast<const irtkMultiLevelTransformation *>(image->Transformation())) ||
  (affd = dynamic_cast<const irtkFreeFormTransformation   *>(image->Transformation()));

  const int nlevels = (mffd ? mffd->NumberOfLevels() : (affd ? 1 : 0));
  if (nlevels == 0) return; // Skip if transformation is not a FFD

  IRTK_START_TIMING();

  for (int lvl = 0; lvl < nlevels; ++lvl) {
    if (mffd) {
      if (!mffd->LocalTransformationIsActive(lvl)) continue;
      affd = mffd->GetLocalTransformation(lvl);
    }

    // Range of control point indices
    blocked_range<int> cps(0, affd->NumberOfCPs());

    // Determine maximum norm of control point gradients
    MaxNodeBasedSimilarityGradient maximum(affd, gradient);
    parallel_reduce(cps, maximum);

    // Sigma value used to suppress noise
    const double sigma = _NodeBasedPreconditioning * maximum.Norm();

    // Normalize control point gradients to be possibly similar
    NormalizeNodeBasedSimilarityGradient norm(affd, gradient, sigma);
    parallel_for(cps, norm);

    // Gradient w.r.t parameters of next active level
    gradient += affd->NumberOfDOFs();
  }

  IRTK_DEBUG_TIMING(2, "normalization of (dis-)similarity gradient");
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::EvaluateGradient(irtkRegisteredImage  *image,
                                           GradientImageType   *&np_gradient,
                                           double               *gradient,
                                           double step, double weight)
{
  if (!image->Transformation()) return;
  const int ndofs = image->Transformation()->NumberOfDOFs();
  // Directly update gradient if no node-based preconditioning is used
  GradientType * const tmp_gradient = _Gradient;
  if (_NodeBasedPreconditioning <= .0) _Gradient = gradient;
  // Compute parametric gradient w.r.t. given transformed image
  if (_Gradient != gradient) {
    memset(_Gradient, 0, ndofs * sizeof(double));
  }
  // Compute analytic gradient if implemented by subclass
  if (!_UseApproximateGradient) {
    if (!np_gradient) {
      np_gradient = new GradientImageType(_Domain, 3);
    }
    if (this->NonParametricGradient(image, np_gradient)) {
      // Normalize voxel-based gradient
      if (_VoxelWisePreconditioning > .0) {
        this->NormalizeGradient(np_gradient);
      }
      this->ParametricGradient(image, np_gradient, _Gradient, weight);
    // ...otherwise, use finite differences approximation
    } else {
      _UseApproximateGradient = true;
    }
  }
  // If no analytic gradient computation implemented or approximation chosen,
  // approximate the image similarity measure gradient using finite differences
  if (_UseApproximateGradient) {
    Delete(np_gradient);
    this->ApproximateGradient(image, _Gradient, step, weight);
  }
  // Normalize node-based gradient
  if (_NodeBasedPreconditioning > .0) {
    this->NormalizeGradient(image, _Gradient);
    for (int dof = 0; dof < ndofs; ++dof) {
      // Note: Normalization cancels out weight factor!
      gradient[dof] += weight * _Gradient[dof];
    }
  }
  // Restore gradient pointer
  _Gradient = tmp_gradient;
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::EvaluateGradient(double *gradient, double step, double weight)
{
  // Compute parametric gradient w.r.t target transformation
  this->EvaluateGradient(_Target, _GradientWrtTarget, gradient, step, weight);
  // If target and source are transformed by different transformations,
  // the gradient vector contains first the derivative values w.r.t the
  // parameters of the target transformation followed by those computed
  // w.r.t the parameters of the source transformation. Otherwise, if
  // both images are transformed by the same transformation, i.e., a
  // velocity based transformation integrated half way in both directions,
  // the derivative values are summed up instead.
  if (_Target->Transformation() && _Source->Transformation()) {
    if (!HaveSameDOFs(_Target->Transformation(), _Source->Transformation())) {
      gradient += _Target->Transformation()->NumberOfDOFs();
    }
  }
  // Compute parametric gradient w.r.t source transformation
  this->EvaluateGradient(_Source, _GradientWrtSource, gradient, step, weight);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageSimilarity::Print(irtkIndent indent) const
{
  irtkEnergyTerm::Print(indent);
  cout << "Target image:";
  if (_Target) {
    cout << endl;
    _Target->Print(indent + 1);
  } else {
    cout << " Missing!" << endl;
  }
  cout << "Source image:";
  if (_Source) {
    cout << endl;
    _Source->Print(indent + 1);
  } else {
    cout << " Missing!" << endl;
  }
  cout << "Image domain:" << endl;
  _Domain.Print(indent + 1);
  cout << "Overlap mask:";
  if (_Mask) {
    cout << endl;
    _Mask->Print(indent + 1);
  } else {
    cout << " None" << endl;
  }
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  if (_Target->Transformation() || all) {
    snprintf(fname, sz, "%starget%s.nii.gz", prefix, suffix);
    _Target->GetFrame(0).Write(fname);
    if (_Target->T() > 4) {
      snprintf(fname, sz, "%starget_gradient%s.nii.gz", prefix, suffix);
      _Target->GetFrame(1, 3).Write(fname);
    }
  }
  if (_Source->Transformation() || all) {
    snprintf(fname, sz, "%ssource%s.nii.gz", prefix, suffix);
    _Source->GetFrame(0).Write(fname);
    if (_Source->T() > 4) {
      snprintf(fname, sz, "%ssource_gradient%s.nii.gz", prefix, suffix);
      _Source->GetFrame(1, 3).Write(fname);
    }
  }
}

// -----------------------------------------------------------------------------
void irtkImageSimilarity::WriteGradient(const char *p, const char *suffix) const
{
  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  // Image derivatives needed for gradient computation
  if (_Target->T() >= 4) {
    snprintf(fname, sz, "%starget_gradient%s.nii.gz", prefix, suffix);
    _Target->GetFrame(1, 3).Write(fname);
  }
  if (_Target->T() > 4) {
    snprintf(fname, sz, "%starget_hessian%s.nii.gz", prefix, suffix);
    _Target->GetFrame(4, _Target->T()-1).Write(fname);
  }

  if (_Source->T() >= 4) {
    snprintf(fname, sz, "%ssource_gradient%s.nii.gz", prefix, suffix);
    _Source->GetFrame(1, 3).Write(fname);
  }
  if (_Source->T() > 4) {
    snprintf(fname, sz, "%ssource_hessian%s.nii.gz", prefix, suffix);
    _Source->GetFrame(4, _Source->T()-1).Write(fname);
  }

  // Computed non-parametric similarity gradient
  if (_GradientWrtTarget) {
    snprintf(fname, sz, "%sgradient_wrt_target%s.nii.gz", prefix, suffix);
    _GradientWrtTarget->Write(fname);
  }
  if (_GradientWrtSource) {
    snprintf(fname, sz, "%sgradient_wrt_source%s.nii.gz", prefix, suffix);
    _GradientWrtSource->Write(fname);
  }
}
