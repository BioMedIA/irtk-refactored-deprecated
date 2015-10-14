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

#include <irtkTransformation.h>

void PrintHelp(const char *name)
{
  cout << "Usage: " << name << " <dofin> <dofout> [options]" << endl;
  PrintStandardOptions(cout);
}

int main(int argc, char **argv)
{
  // Parse arguments
  REQUIRES_POSARGS(2);

  for (ALL_OPTIONS) {
    HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  // Read transformation
  auto_ptr<irtkAffineTransformation> transformation(new irtkAffineTransformation);
  transformation->Read(POSARG(1));
  if (verbose) transformation->Print();

  // Invert transformation
  if (verbose) cout << "Inverting transformation ...";
  transformation->Invert();
  if (verbose) cout << " done" << endl;

  // Write transform
  transformation->Write(POSARG(2));
  if (verbose) transformation->Print();

  return 0;
}
