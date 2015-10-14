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

#include <irtkDefines.h> // WINDOWS
#include <irtkCxxLib.h>
#include <irtkPath.h>


// -----------------------------------------------------------------------------
#if WINDOWS
const char PATHSEP = '\\';
#else
const char PATHSEP = '/';
#endif

// -----------------------------------------------------------------------------
string Extension(const char *path, ExtensionMode mode)
{
  // Get base name without directory path
  string s = BaseName(path);
  if (s.empty() || mode == EXT_None) return "";
  // Split file path, note that first part is the file name itself and that
  // an empty initial part is caused by the leading '.' of hidden files on Unix
  vector<string> parts = Split(path, ".");
#if WINDOWS
  const bool hidden = false;
#else
  const bool hidden = (s[0] == '.');
#endif
  parts.erase(parts.begin(), parts.begin() + (hidden ? 2 : 1));
  // If no parts left, return empty string
  if (parts.empty()) return "";
  // Concatenate requested file extension parts only
  switch (mode) {
    case EXT_Default:
    case EXT_LastWithGz:
      s  = '.';
      s += parts.back();
      transform(s.begin(), s.end(), s.begin(), ::tolower);
      if (s == ".gz" && parts.size() > 1) s = parts[parts.size()-2] + '.' + s;
      break;
    case EXT_LastWithoutGz:
      s  = '.';
      s += parts.back();
      transform(s.begin(), s.end(), s.begin(), ::tolower);
      if (s == ".gz" && parts.size() > 1) s = '.' + parts[parts.size()-2];
      break;
    case EXT_Last:
      s  = '.';
      s += parts.back();
      break;
    case EXT_All:
      s.clear();
      for (size_t i = 0; i < parts.size(); ++i) {
        s += '.';
        s += parts[i];
      }
      break;
    case EXT_None:
      // Dealt with before
      return "";
  }
  transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}

// -----------------------------------------------------------------------------
string Extension(const string &path, ExtensionMode mode)
{
  return Extension(path.c_str(), mode);
}

// -----------------------------------------------------------------------------
string Directory(const char *path)
{
  string dir(path);
  size_t n = dir.find_last_of(PATHSEP);
  if (n == string::npos) dir.clear();
  else                   dir.resize(n);
  return dir;
}

// -----------------------------------------------------------------------------
string Directory(const string &path)
{
  return Directory(path.c_str());
}

// -----------------------------------------------------------------------------
string BaseName(const char *path)
{
  const char *name = path;
  while (name[0] != '\0') name++;
  while (name != path && name[0] != PATHSEP) name--;
  if (name[0] == PATHSEP) name++;
  return name;
}

// -----------------------------------------------------------------------------
string BaseName(const string &path)
{
  return BaseName(path.c_str());
}

// -----------------------------------------------------------------------------
string FileName(const char *path, ExtensionMode mode)
{
  string name = BaseName(path);
  string ext  = Extension(path, mode);
  return name.substr(0, name.length() - ext.length());
}

// -----------------------------------------------------------------------------
string FileName(const string &path, ExtensionMode mode)
{
  return FileName(path.c_str(), mode);
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
// dirname copied from glibc-2.1.1 with slight changes
char *dirname2(char *path)
{
  static const char dot[] = ".";
  char *last_slash;

  /* Find last '/'.  */
  last_slash = path != NULL ? strrchr(path, PATHSEP) : NULL;

  if (last_slash == path)
  /* The last slash is the first character in the string.  We have to
   return "/".  */
    ++last_slash;
  else if (last_slash != NULL && last_slash[1] == '\0')
  /* The '/' is the last character, we have to look further.  */
    last_slash = (char *)memchr(path, last_slash - path, PATHSEP);

  if (last_slash != NULL)
  /* Terminate the path.  */
    last_slash[0] = '\0';
  else
  /* This assignment is ill-designed but the XPG specs require to
   return a string containing "." in any case no directory part is
   found and so a static and constant string is required.  */
    path = (char *) dot;

  return path;
}

// -----------------------------------------------------------------------------
// basename copied from glibc-2.1.1 with slight changes
char *basename2(char *filename)
{
  char *p = strrchr(filename, PATHSEP);
  return p ? p + 1 : filename;
}
