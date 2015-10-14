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

#ifndef _IRTKPATH_H
#define _IRTKPATH_H

#include <string>


// =============================================================================
// Global constants
// =============================================================================

/// Path separating character
extern const char PATHSEP;

// =============================================================================
// Split file path
// =============================================================================

/// Enumeration of file path extension retrival modes
enum ExtensionMode
{
  EXT_Default,       ///< Default extension mode
  EXT_None,          ///< No part of the file is considered to be an extension
  EXT_Last,          ///< Last file extension only
  EXT_LastWithGz,    ///< Last file extension possibly plus ".gz"
  EXT_LastWithoutGz, ///< Last file extension with possibly ".gz" removed first
  EXT_All            ///< All file extensions, i.e., everything after first dot (besides leading dot of hidden files on Unix)
};

/// Get file name extension in lower case incl. leading dot ('.')
std::string Extension(const char *, ExtensionMode = EXT_Default);

/// Get file name extension in lower case incl. leading dot ('.')
std::string Extension(const std::string &, ExtensionMode = EXT_Default);

/// Get directory part of file path
std::string Directory(const char *);

/// Get directory part of file path
std::string Directory(const std::string &);

/// Get file name of file path incl. file extension
std::string BaseName(const char *);

/// Get file name of file path incl. file extension
std::string BaseName(const std::string &);

/// Get file name of file path excl. file extension
std::string FileName(const char *, ExtensionMode = EXT_Default);

/// Get file name of file path excl. file extension
std::string FileName(const std::string &, ExtensionMode = EXT_Default);

// =============================================================================
// Deprecated
// =============================================================================

/// Get directory component of file path
///
/// @attention Modifies input string by replacing the last PATHSEP by a
///            terminating null character. To get both file path parts,
///            call basename2 prior to dirname2!
///
/// @deprecated Use Directory instead.
char *dirname2(char *);

/// Get file name component of file path including extension
///
/// @deprecated Use BaseName instead.
char *basename2(char *);


#endif
