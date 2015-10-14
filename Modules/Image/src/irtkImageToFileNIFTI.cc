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

#ifdef HAS_NIFTI

#include <irtkImage.h>

#include <irtkImageToFile.h>

#include <irtkNIFTI.h>

// -----------------------------------------------------------------------------
irtkImageToFileNIFTI::irtkImageToFileNIFTI()
:
  _headername(NULL)
{
}

// -----------------------------------------------------------------------------
irtkImageToFileNIFTI::~irtkImageToFileNIFTI()
{
	free(_headername);
}

// -----------------------------------------------------------------------------
// Change memory layout from xyz...xyz...xyz... to xxx...yyy...zzz......
template <class VectorType>
void *ReformatImageData(const void *input, int nx, int ny, int nz, int nt, int nu)
{
  typedef typename voxel_info<VectorType>::ScalarType ScalarType;
  void *output = malloc(nx * ny * nz * nt * nu * sizeof(ScalarType));
  const ScalarType *in  = reinterpret_cast<const ScalarType *>(input);
  ScalarType       *out = reinterpret_cast<      ScalarType *>(output);
  const int offset = nx * ny * nz * nt;
  for (int l = 0; l < nt; ++l) {
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          ScalarType *tmp = out;
          for (int m = 0; m < nu; ++m) {
            *tmp = (*in++);
            tmp += offset;
          }
          ++out;
        }
      }
    }
  }
  return output;
}

// -----------------------------------------------------------------------------
void irtkImageToFileNIFTI::SetOutput(const char *name)
{
  // TODO: Hack to force the use of .nii(.gz) instead of .hdr(.gz)/.img(.gz)?!?
  //       Find out what it is good for and remove this nasty hack if possible.
	free(_headername);
  _headername = strdup(name);
  int length = strlen(_headername);
  if (length > 3) {
    if (strcmp(_headername+length-3, ".GZ") == 0) {
      _headername[length-2] = 'g';
      _headername[length-1] = 'z';
    }
    if (strcmp(_headername+length-3, ".gz") == 0) {
      length -= 3;
    }
  }
  if (length > 3) {
    _headername[length-3] = 'n';
    _headername[length-2] = 'i';
    _headername[length-1] = 'i';
  }
	irtkImageToFile::SetOutput(_headername);
}

// -----------------------------------------------------------------------------
void irtkImageToFileNIFTI::Initialize()
{
  // Get image info
  const int nx = _input->X();
  const int ny = _input->Y();
  const int nz = _input->Z();
  const int nt = _input->T();

  double xsize, ysize, zsize, tsize;
  _input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);

  irtkMatrix qmat = _input->GetImageToWorldMatrix();
  irtkMatrix *smat = NULL;
  if (!_input->GetAffineMatrix().IsIdentity()) {
    smat = new irtkMatrix(qmat);
    qmat = _input->GetAffineMatrix().Inverse() * qmat;
  }
  const double torigin = _input->ImageToTime(0);

  // Init header
  void * const data = _input->GetDataPointer();
  switch (_input->GetDataType()) {
    case IRTK_VOXEL_UNSIGNED_CHAR: {
      _hdr.Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_UINT8, qmat, smat, torigin, data);
      break;
    }
    case IRTK_VOXEL_SHORT: {
      _hdr.Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_INT16, qmat, smat, torigin, data);
      break;
    }
    case IRTK_VOXEL_UNSIGNED_SHORT: {
      _hdr.Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_UINT16, qmat, smat, torigin, data);
      break;
    }
    case IRTK_VOXEL_INT: {
      _hdr.Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_INT32, qmat, smat, torigin, data);
      break;
    }
    case IRTK_VOXEL_UNSIGNED_INT: {
      _hdr.Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_UINT32, qmat, smat, torigin, data);
      break;
    }
    case IRTK_VOXEL_FLOAT: {
      _hdr.Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_FLOAT32, qmat, smat, torigin, data);
      break;
    }
//    case IRTK_VOXEL_FLOAT2: {
//      _hdr.Initialize(nx, ny, nz, nt, 2, xsize, ysize, zsize, tsize,
//                      NIFTI_TYPE_FLOAT32, qmat, smat, torigin,
//                      ReformatImageData<irtkFloat2>(data, nx, ny, nz, nt, 2));
//      break;
//    }
    case IRTK_VOXEL_FLOAT3: {
      _hdr.Initialize(nx, ny, nz, nt, 3, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_FLOAT32, qmat, smat, torigin,
                      ReformatImageData<irtkFloat3>(data, nx, ny, nz, nt, 3));
      break;
    }
    case IRTK_VOXEL_FLOAT4: {
      _hdr.Initialize(nx, ny, nz, nt, 4, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_FLOAT32, qmat, smat, torigin,
                      ReformatImageData<irtkFloat4>(data, nx, ny, nz, nt, 4));
      break;
    }
    case IRTK_VOXEL_DOUBLE: {
      _hdr.Initialize(nx, ny, nz, nt, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_FLOAT64, qmat, smat, torigin, data);
      break;
    }
//    case IRTK_VOXEL_DOUBLE2: {
//      _hdr.Initialize(nx, ny, nz, nt, 2, xsize, ysize, zsize, tsize,
//                      NIFTI_TYPE_FLOAT64, qmat, smat, torigin,
//                      ReformatImageData<irtkDouble2>(data, nx, ny, nz, nt, 2));
//      break;
//    }
    case IRTK_VOXEL_DOUBLE3: {
      _hdr.Initialize(nx, ny, nz, nt, 3, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_FLOAT64, qmat, smat, torigin,
                      ReformatImageData<irtkDouble3>(data, nx, ny, nz, nt, 3));
      break;
    }
    case IRTK_VOXEL_DOUBLE4: {
      _hdr.Initialize(nx, ny, nz, nt, 4, xsize, ysize, zsize, tsize,
                      NIFTI_TYPE_FLOAT64, qmat, smat, torigin,
                      ReformatImageData<irtkDouble4>(data, nx, ny, nz, nt, 4));
      break;
    }
    default:
      cerr << "irtkImageToFileNIFTI::Run(): Unknown voxel type" << endl;
      exit(1);
  }

  delete smat;
}

// -----------------------------------------------------------------------------
void irtkImageToFileNIFTI::Finalize()
{
  if (_hdr.nim->data != _input->GetScalarPointer()) free(_hdr.nim->data);
  _hdr.nim->data = NULL;
}

// -----------------------------------------------------------------------------
void irtkImageToFileNIFTI::Run()
{
	// Initialize filter
	this->Initialize();

	// Set filename in _hdr
	struct nifti_1_header nhdr = nifti_convert_nim2nhdr(_hdr.nim);

  void * const data = _hdr.nim->data; // Keep copy of data pointer
  _hdr.nim->data    = NULL;           // Free nifti_image, but not data
  nifti_image_free(_hdr.nim);

	_hdr.nim               = nifti_convert_nhdr2nim(nhdr, _output); // This sets fname and iname
	_hdr.nim->iname_offset = 352;       // Some nifti versions lose this on the way!
  _hdr.nim->data         = data;      // Restore data pointer

	// Write hdr and data
	nifti_image_write(_hdr.nim);

	// Finalize filter
	this->Finalize();
}


#endif
