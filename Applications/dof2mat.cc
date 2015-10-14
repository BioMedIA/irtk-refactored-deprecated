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

#include <irtkTransformation.h>

char *matout_name = NULL;

void usage()
{
  cerr << "Usage: dof2mat [doffile] [-matout matrixfile] [-invert]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  irtkMatrix matrix;
  irtkAffineTransformation transformation;

  if (argc < 2) {
    usage();
  }

  // Read transformation
  transformation.irtkTransformation::Read(argv[1]);
  argc--;
  argv++;

  // Convert to matrix
  matrix = transformation.GetMatrix();

// Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-matout") == 0)) {
      argc--;
      argv++;
      matout_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      matrix.Invert();
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Print matrix
  matrix.Print();

  // Write matrix if necessary
  if (matout_name != NULL) {
    matrix.Write(matout_name);
  }
}
