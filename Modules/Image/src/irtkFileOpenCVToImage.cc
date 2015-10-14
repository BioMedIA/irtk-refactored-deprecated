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

#ifdef HAS_OPENCV

#include <irtkImageToOpenCv.h>

int irtkFileOpenCVToImage::CheckHeader(const char *filename)
{
    IplImage *pimage = NULL;
    pimage= cvLoadImage(filename,0);
    if(pimage != NULL){
        cvReleaseImage(&pimage);
        return 1;
    }else{
        return 0;
    }
}

void irtkFileOpenCVToImage::SetInput(const char *filename){
    _pimage = NULL;
    _pimage= cvLoadImage(filename,0);

    // Read header
    this->ReadHeader();
}

irtkImage * irtkFileOpenCVToImage::GetOutput(){
    irtkGreyImage *output;
    irtkImageToOpenCv<irtkGreyPixel> itocg;
    itocg.SetOutput(_pimage);
    itocg.Invert();
    output = itocg.GetInput();
    return (irtkImage*)output;
}

void irtkFileOpenCVToImage::ReadHeader()
{
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
}

#endif

