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

#include <irtkGenericRegistrationFilter.h>
#include <irtkGenericRegistrationDebugger.h>

#include <irtkTransformation.h>
#include <irtkImageSimilarity.h>

#ifdef HAS_VTK
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkShortArray.h>
#include <vtkFloatArray.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <irtkPolyDataUtils.h>
using namespace irtk::polydata;
#endif

// -----------------------------------------------------------------------------
template <class TReal>
void WriteGradient(const char *fname, irtkFreeFormTransformation *ffd, int l, const TReal *g)
{
  irtkGenericImage<TReal> gradient(ffd->Attributes(), 3);
  int xdof, ydof, zdof;
  for (int k = 0; k < ffd->Z(); ++k) {
    for (int j = 0; j < ffd->Y(); ++j) {
      for (int i = 0; i < ffd->X(); ++i) {
        ffd->IndexToDOFs(ffd->LatticeToIndex(i, j, k, l), xdof, ydof, zdof);
        gradient(i, j, k, 0) = g[xdof];
        gradient(i, j, k, 1) = g[ydof];
        gradient(i, j, k, 2) = g[zdof];
      }
    }
  }
  gradient.Write(fname);
}

// -----------------------------------------------------------------------------
static void WriteAsVTKDataSet(const char *fname, irtkFreeFormTransformation *ffd, const double *g = NULL)
{
#ifdef HAS_VTK
  vtkSmartPointer<vtkPoints>         pos  = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkShortArray>     stat = vtkSmartPointer<vtkShortArray>::New();
  vtkSmartPointer<vtkFloatArray>     coef = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkFloatArray>     disp = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkFloatArray>     grad = (g ? vtkSmartPointer<vtkFloatArray>::New() : NULL);
  vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

  pos->SetNumberOfPoints(ffd->NumberOfCPs());

  stat->SetName("status");
  stat->SetNumberOfComponents(1);
  stat->SetNumberOfTuples(ffd->NumberOfCPs());

  coef->SetName("coefficient");
  coef->SetNumberOfComponents(3);
  coef->SetNumberOfTuples(ffd->NumberOfCPs());

  disp->SetName("displacement");
  disp->SetNumberOfComponents(3);
  disp->SetNumberOfTuples(ffd->NumberOfCPs());

  if (grad) {
    grad->SetName("gradient");
    grad->SetNumberOfComponents(3);
    grad->SetNumberOfTuples(ffd->NumberOfCPs());
  }

  int    i, j, k;
  double x1, y1, z1, x2, y2, z2;
  for (int cp = 0; cp < ffd->NumberOfCPs(); ++cp) {
    ffd->IndexToLattice(cp, i, j, k);
    x1 = i, y1 = j, z1 = k;
    ffd->LatticeToWorld(x1, y1, z1);
    pos->SetPoint(cp, x1, y1, z1);
    x2 = x1, y2 = y1, z2 = z1;
    ffd->Transform(x2, y2, z2);
    disp->SetTuple3(cp, x2 - x1, y2 - y1, z2 - z1);
    stat->SetTuple1(cp, ffd->IsActive(cp));
    ffd->Get(i, j, k, x2, y2, z2);
    coef->SetTuple3(cp, x2, y2, z2);
    if (grad) {
      ffd->IndexToDOFs(cp, i, j, k);
      grad->SetTuple3(cp, g[i], g[j], g[k]);
    }
  }

  grid->SetDimensions(ffd->X(), ffd->Y(), ffd->Z());
  grid->SetPoints(pos);
  grid->GetPointData()->SetScalars(stat);
  grid->GetPointData()->SetVectors(coef);
  grid->GetPointData()->AddArray(disp);
  if (grad) grid->GetPointData()->AddArray(grad);

  vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(fname);
  writer->SetCompressorTypeToZLib();
  SetVTKInput(writer, grid);
  writer->Update();
#endif // defined(HAS_VTK)
}

// -----------------------------------------------------------------------------
void WriteTransformation(const char                    *fname,
                         irtkHomogeneousTransformation *lin,
                         const irtkPoint               &target_offset,
                         const irtkPoint               &source_offset)
{
  const irtkMatrix mat = lin->GetMatrix();
  irtkMatrix pre (4, 4);
  irtkMatrix post(4, 4);
  pre .Ident();
  post.Ident();
  pre (0, 3) = - target_offset._x;
  pre (1, 3) = - target_offset._y;
  pre (2, 3) = - target_offset._z;
  post(0, 3) = + source_offset._x;
  post(1, 3) = + source_offset._y;
  post(2, 3) = + source_offset._z;
  lin->PutMatrix(post * mat * pre);
  lin->Write(fname);
  lin->PutMatrix(mat);
}

// -----------------------------------------------------------------------------
irtkGenericRegistrationDebugger::irtkGenericRegistrationDebugger(const char *prefix)
:
  irtkObserver(),
  _Prefix      (prefix),
  _LevelPrefix (true),
  _Registration(NULL)
{
}

// -----------------------------------------------------------------------------
irtkGenericRegistrationDebugger::~irtkGenericRegistrationDebugger()
{
}

// -----------------------------------------------------------------------------
void irtkGenericRegistrationDebugger::HandleEvent(irtkObservable *obj, irtkEvent event, const void *data)
{
  irtkGenericRegistrationFilter * const r = _Registration;

  const int sz = 256;
  char prefix[sz];
  char suffix[sz];
  char fname [sz];

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Initialize/update debugger, set common file name prefix/suffix
  switch (event) {

    // -------------------------------------------------------------------------
    // Attach/detach debugger
    case RegisteredEvent:
      _Registration = dynamic_cast<irtkGenericRegistrationFilter *>(obj);
      if (!_Registration) {
        cerr << "irtkGenericRegistrationDebugger::HandleEvent: Cannot attach debugger to object which is not of type irtkGenericRegistrationFilter" << endl;
        exit(1);
      }
      _Iteration = 0;
      break;
    case UnregisteredEvent:
      _Registration = NULL;
      break;

    // -------------------------------------------------------------------------
    // Start/end
    case StartEvent:
      _Level = reinterpret_cast<const irtkIteration *>(data)->Count() + 1; // Iter()
      if (_LevelPrefix) _Iteration = 0;
      // Get pointers to similarity terms
      _Similarity.clear();
      for (int i = 0; i < r->_Energy.NumberOfTerms(); ++i) {
        irtkImageSimilarity *similarity = dynamic_cast<irtkImageSimilarity *>(r->_Energy.Term(i));
        if (similarity) _Similarity.push_back(similarity);
      }
      // Do not add a break statement here!
    case EndEvent:
      snprintf(prefix, sz, "%slevel_%d_", _Prefix.c_str(), _Level);
      suffix[0] = '\0';
      break;

    // -------------------------------------------------------------------------
    // Iteration
    case IterationEvent:
      ++_Iteration;
      return; // No data to write yet
    case LineSearchIterationStartEvent:
      _LineIteration = reinterpret_cast<const irtkIteration *>(data)->Iter();
      return; // No data to write yet

    case LineSearchStartEvent:
      if (_LevelPrefix) {
        snprintf(prefix, sz, "%slevel_%d_",         _Prefix.c_str(), _Level);
        snprintf(suffix, sz, "_%03d",               _Iteration);
      } else {
        strcpy(prefix, _Prefix.c_str());
        snprintf(suffix, sz, "_%03d",               _Iteration);
      }
      break;
    case AcceptedStepEvent:
      if (_LevelPrefix) {
        snprintf(prefix, sz, "%slevel_%d_",         _Prefix.c_str(), _Level);
        snprintf(suffix, sz, "_%03d_%03d_accepted", _Iteration, _LineIteration);
      } else {
        strcpy(prefix, _Prefix.c_str());
        snprintf(suffix, sz, "_%03d_%03d_accepted", _Iteration, _LineIteration);
      }
      break;
    case RejectedStepEvent:
      if (_LevelPrefix) {
        snprintf(prefix, sz, "%slevel_%d_",         _Prefix.c_str(), _Level);
        snprintf(suffix, sz, "_%03d_%03d_rejected", _Iteration, _LineIteration);
      } else {
        strcpy(prefix, _Prefix.c_str());
        snprintf(suffix, sz, "_%03d_%03d_rejected", _Iteration, _LineIteration);
      }
      break;

    // -------------------------------------------------------------------------
    // Ignored event
    default: return;
  }

  irtkMultiLevelTransformation  *mffd = NULL;
  irtkFreeFormTransformation    *ffd  = NULL;
  irtkHomogeneousTransformation *lin  = NULL;

  if (r) {
    (mffd = dynamic_cast<irtkMultiLevelTransformation  *>(r->_Transformation)) ||
    (ffd  = dynamic_cast<irtkFreeFormTransformation    *>(r->_Transformation)) ||
    (lin  = dynamic_cast<irtkHomogeneousTransformation *>(r->_Transformation));
  }
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (mffd->LocalTransformationIsActive(l)) {
        if (ffd) {
          ffd = NULL;
          break;
        }
        ffd = mffd->GetLocalTransformation(l);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Write debug information
  switch (event) {

    // -------------------------------------------------------------------------
    // Write initial state
    case StartEvent: {

      // Write input images and their derivatives
      for (size_t i = 0; i < r->_Image[r->_CurrentLevel].size(); ++i) {
        snprintf(fname, sz, "%simage_%02lu.nii.gz", prefix, i+1);
        r->_Image[r->_CurrentLevel][i].Write(fname);
        if (debug >= 2) {
          irtkBaseImage *gradient = NULL;
          irtkBaseImage *hessian  = NULL;
          for (size_t j = 0; j < _Similarity.size(); ++j) {
            if (_Similarity[j]->Target()->InputImage() == &r->_Image[r->_CurrentLevel][i]) {
              gradient = _Similarity[j]->Target()->InputGradient();
              hessian  = _Similarity[j]->Target()->InputHessian();
              break;
            }
            if (_Similarity[j]->Source()->InputImage() == &r->_Image[r->_CurrentLevel][i]) {
              gradient = _Similarity[j]->Source()->InputGradient();
              hessian  = _Similarity[j]->Source()->InputHessian();
              break;
            }
          }
          if (gradient) {
            snprintf(fname, sz, "%simage_%02lu_gradient.nii.gz", prefix, i+1);
            gradient->Write(fname);
          }
          if (hessian) {
            snprintf(fname, sz, "%simage_%02lu_hessian.nii.gz", prefix, i+1);
            hessian->Write(fname);
          }
        }
      }

#ifdef HAS_VTK
      // Write input point set
      for (size_t i = 0; i < r->_PointSet[r->_CurrentLevel].size(); ++i) {
        vtkPointSet *pointset = r->_PointSet[r->_CurrentLevel][i];
        snprintf(fname, sz, "%spointset_%02lu%s", prefix, i+1, DefaultExtension(pointset));
        WritePointSet(fname, pointset);
      }
#endif

    } break;

    // -------------------------------------------------------------------------
    // Write intermediate results after each gradient step
    case LineSearchStartEvent: {

      // Energy gradient vector
      const double * const gradient = reinterpret_cast<const irtkLineSearchStep *>(data)->_Direction;

      // Write input and other debug output of energy terms
      r->_Energy.WriteDataSets(prefix, suffix, _Iteration == 1);

      if (debug >= 3) {

        // Write non-parametric gradient of data fidelity terms
        r->_Energy.WriteGradient(prefix, suffix);

        // Write energy gradient(s) w.r.t control points
        if (ffd) {
          if (ffd->T() > 1) {
            for (int l = 0; l < ffd->T(); ++l) {
              snprintf(fname, sz, "%senergy_gradient_t%02d%s.nii.gz", prefix, l+1, suffix);
              WriteGradient(fname, ffd, l, gradient);
            }
          } else {
            snprintf(fname, sz, "%senergy_gradient%s.nii.gz", prefix, suffix);
            WriteGradient(fname, ffd, 0, gradient);
          }
        } else if (mffd) {
          const double *g = gradient;
          for (int i = 0; i < mffd->NumberOfLevels(); ++i) {
            if (!mffd->LocalTransformationIsActive(i)) continue;
            ffd = mffd->GetLocalTransformation(i);
            if (ffd->T() > 1) {
              for (int l = 0; l < ffd->T(); ++l) {
                snprintf(fname, sz, "%senergy_gradient_wrt_ffd_%d_t%02d%s.nii.gz", prefix, i+1, l+1, suffix);
                WriteGradient(fname, ffd, l, gradient);
              }
            } else {
              snprintf(fname, sz, "%senergy_gradient_wrt_ffd_%d_%s.nii.gz", prefix, i+1, suffix);
              WriteGradient(fname, ffd, 0, g);
            }
            g += ffd->NumberOfDOFs();
          }
          ffd = NULL;
        } else if (lin) {
          snprintf(fname, sz, "%senergy_gradient%s.txt", prefix, suffix);
          ofstream of(fname);
          for (int dof = 0; dof < r->_Energy.NumberOfDOFs(); ++dof) {
            of << gradient[dof] << "\n";
          }
          of.close();
        }

      }

      // Write current transformation estimate
      snprintf(fname, sz, "%stransformation%s.dof.gz", prefix, suffix);
      if (lin) {
        WriteTransformation(fname, lin, r->_TargetOffset, r->_SourceOffset);
      } else {
        r->_Transformation->Write(fname);
        if (ffd && r->_Input.empty() && r->NumberOfPointSets() > 0 && debug >= 4) {
          snprintf(fname, sz, "%stransformation%s.vtp", prefix, suffix);
          WriteAsVTKDataSet(fname, ffd, gradient);
        }
      }
    } break;

    case AcceptedStepEvent:
    case RejectedStepEvent: {
      if (debug >= 5) {

        // Write updated input of data fidelity terms
        r->_Energy.WriteDataSets(prefix, suffix, false);

        // Write current transformation estimate
        snprintf(fname, sz, "%stransformation%s.dof.gz", prefix, suffix);
        if (lin) WriteTransformation(fname, lin, r->_TargetOffset, r->_SourceOffset);
        else     r->_Transformation->Write(fname);

      }
    } break;

    // -------------------------------------------------------------------------
    // Unhandled event
    default: break;
  }
}
