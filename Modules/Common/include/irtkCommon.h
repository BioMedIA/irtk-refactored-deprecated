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

#ifndef _IRTKCOMMON_H
#define _IRTKCOMMON_H


#include <irtkDefines.h>
#if WINDOWS
#  include <irtkWindows.h>
#endif
#include <irtkCxxLib.h>
#include <irtkZLib.h>
#include <irtkObject.h>
#include <irtkVersion.h>
#include <irtkObserver.h>
#include <irtkObservable.h>
#include <irtkCifstream.h>
#include <irtkCofstream.h>
#include <irtkIndent.h>
#include <irtkAllocate.h>
#include <irtkDeallocate.h>
#include <irtkException.h>
#include <irtkOptions.h>
#include <irtkParallel.h>
#include <irtkPath.h>
#include <irtkProfiling.h>
#include <irtkFloat.h>
#include <irtkMath.h>
#include <irtkEnums.h>
#include <irtkTerminal.h>
#include <irtkTestProd.h>

#include "FastDelegate.h"

// TODO: Consider moving these to the PolyData package
//#ifdef HAS_VTK
//#include <irtkPolyDataUtils.h>
//using namespace irtk::polydata;
//#endif


/// Median of weighted sequence
double weightedmedian(int, double,double, float*, float*);


#endif
