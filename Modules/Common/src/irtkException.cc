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

#include <irtkCommon.h>

std::ostream& operator<<(std::ostream& os, const irtkException& ex)
{
  os << "File name: ";
  if (ex.GetFileName() != "")
    os << ex.GetFileName();
  else
    os << "UNKNOWN";
  os << std::endl;

  os << "Line: ";
  if (ex.GetLine() != 0)
    os << ex.GetLine();
  else
    os << "UNKNOWN";
  os << std::endl;

  os << "Message: " << ex.GetMessage();

  return os;
}
