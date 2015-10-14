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

#ifndef _IRTKOPTIONS_H
#define _IRTKOPTIONS_H

#include <iostream>


/*
   \file  irtkOptions.h
   \brief Command-line parsing library.

   The macros defined by this module should reduce the number of lines of code
   required to parse simple command-line arguments in the main function of a
   program. Use them for example as follows:
  
   \code
   for (ALL_OPTIONS) {
     if      (OPTION("-steps" )) NumberOfIntegrationSteps = atoi(ARGUMENT);
     else if (OPTION("-iters" )) NumberOfBCHSteps         = atoi(ARGUMENT);
     else if (OPTION("-terms" )) NumberOfBCHTerms         = atoi(ARGUMENT);
     else if (OPTION("-smooth")) SmoothBCHApproximation   = true;
     else HANDLE_COMMON_OR_UNKNOWN_OPTION();
   }
   \endcode
  
   If your program has positional arguments which have to be specified by the
   user before the optional arguments, use the following code instead:
  
   \code
   REQUIRES_POSARGS(2);
   const char *input_file  = POSARG(1);
   const char *output_file = POSARG(2);
   for (ALL_OPTIONS) {
     // same as above
   }
   \endcode
  
   Other uses including options with multiple arguments may look like this:
   \code
   REQUIRES_POSARGS(0);
   DISCARD_PARSED_OPTIONS();
   int    x  = 1, y  = 1, z  = 1, t  = 1;
   double dx = 1, dy = 1, dz = 1, dt = 1;
   for (ALL_OPTIONS) {
     // -dim <x> <y> <z> <t>
     if (OPTION("-dim")) {
       x = atoi(ARGUMENT);
       y = atoi(ARGUMENT);
       z = atoi(ARGUMENT);
       t = atoi(ARGUMENT);
     } else if (OPTION("-pixdim")) {
       dx = atof(ARGUMENT);
       dy = atof(ARGUMENT);
       dz = atof(ARGUMENT);
       dt = atof(ARGUMENT);
     } else HANDLE_COMMON_OPTION();
   }
   // do something, then parse remaining options
   bool foo = false;
   for (ALL_OPTIONS) {
     if (OPTION("-foo")) foo = true;
     else HANDLE_UNKNOWN_OPTION();
   }
   \endcode
 */

// =============================================================================
// Global standard options
// =============================================================================

class irtkVersion;

/// Current software version
extern const irtkVersion current_version;

/// Version to emulate
extern irtkVersion version;

/// Verbosity of output messages
extern int verbose;

/// Debug level, e.g., amount of intermediate data to write to disk
/// This flag can be combined with verbose after parsing the command options.
/// For example, "if (debug) verbose = 100;", or treated separately.
extern int debug;

// =============================================================================
// Standard options
// =============================================================================

/// Print software revision number (or version if not available) only
void PrintRevision(std::ostream &);

/// Print build time stamp as version string
void PrintVersion(std::ostream &, const char* = NULL);

/// Check if given option is a standard option
bool IsStandardOption(const char *);

/// Parse standard option
void ParseStandardOption(int &, int &, char *[]);

//// Print standard options of any IRTK command
void PrintStandardOptions(std::ostream &);

//// Print common options of any IRTK command
void PrintCommonOptions(std::ostream &);

// =============================================================================
// Command helper macros
// =============================================================================

// -----------------------------------------------------------------------------
#define EXECNAME                   const_cast<const char *>(basename2(const_cast<char *>(argv[0])))
#define ARGIDX                     _i
#define OPTIDX                     ARGIDX
#define OPTNAME                    _option
#define DISCARD_PARSED_POSARGS()   _discard_parsed_posargs = true
#define DISCARD_PARSED_OPTIONS()   _discard_parsed_options = true
#define DISCARD_PARSED_ARGUMENTS() _discard_parsed_options = _discard_parsed_posargs = true
#define KEEP_PARSED_POSARGS()      _discard_parsed_posargs = false
#define KEEP_PARSED_OPTIONS()      _discard_parsed_options = false
#define KEEP_PARSED_ARGUMENTS()    _discard_parsed_options = _discard_parsed_posargs = false
#define ALL_ARGUMENTS              int ARGIDX = 1; ARGIDX < argc; ARGIDX++
#define ALL_POSARGS                int ARGIDX = 1; ARGIDX <= _numposarg; ARGIDX++
#define OPTIONAL_POSARGS           int ARGIDX = _posargc + 1; ARGIDX < argc && argv[ARGIDX][0] != '-'; ARGIDX++
#define ARGUMENTS_AFTER(pos)       int ARGIDX = (pos)+1; ARGIDX < argc; ARGIDX++
#define ALL_OPTIONS                ARGUMENTS_AFTER(_numposarg == -1 ? _posargc : _numposarg)
#define POSARG(i)                  _GetPositionalArgument((i), argc, argv)
#define NUM_POSARGS                (_numposarg == -1 ? _GetNumberOfPositionalArguments(argc, argv) : _numposarg)
#define IS_OPTION                  _IsOption(OPTIDX, argc, argv)
#define OPTION(opt)                _IsOption(OPTIDX, argc, argv, opt)
#define ARGUMENT                   _GetOptionArgument(OPTIDX, argc, argv)
#define HAS_ARGUMENT               _IsArgument(OPTIDX, argc, argv)
#define HELP_OPTION                (OPTION("-h") || OPTION("-help") || OPTION("--help"))
#define VERSION_OPTION             (OPTION("-version") || OPTION("--version") || OPTION("-revision"))
#define STANDARD_OPTION            IsStandardOption    (argv[OPTIDX])
#define PARALLEL_OPTION            IsParallelOption    (argv[OPTIDX])
#define PROFILING_OPTION           IsProfilingOption   (argv[OPTIDX])
#define TERMINAL_OPTION            IsTerminalOption    (argv[OPTIDX])
#define PARSE_STANDARD_OPTION()    ParseStandardOption (OPTIDX, argc, argv)
#define PARSE_PARALLEL_OPTION()    ParseParallelOption (OPTIDX, argc, argv)
#define PARSE_PROFILING_OPTION()   ParseProfilingOption(OPTIDX, argc, argv)
#define PARSE_TERMINAL_OPTION()    ParseTerminalOption (OPTIDX, argc, argv)

// -----------------------------------------------------------------------------
#define REQUIRES_POSARGS(n)                                                    \
  do {                                                                         \
    _posargc = n;                                                              \
    for (ALL_ARGUMENTS) {                                                      \
      if      (HELP_OPTION   ) HANDLE_HELP_OPTION();                           \
      else if (VERSION_OPTION) HANDLE_VERSION_OPTION();                        \
    }                                                                          \
    if (argc <= _posargc) {                                                    \
      PrintHelp(EXECNAME);                                                     \
      exit(1);                                                                 \
    }                                                                          \
    _numposarg = NUM_POSARGS;                                                  \
    if (_numposarg < _posargc) _numposarg = _posargc;                          \
  } while (false)

// -----------------------------------------------------------------------------
#define EXPECTS_POSARGS(n)                                                     \
  REQUIRES_POSARGS(n);                                                         \
  do {                                                                         \
    if (NUM_POSARGS > n) {                                                     \
      PrintHelp(EXECNAME);                                                     \
      exit(1);                                                                 \
    }                                                                          \
  } while (false)

// -----------------------------------------------------------------------------
#define HANDLE_HELP_OPTION()                                                   \
  do { PrintHelp(EXECNAME); exit(0); } while (false)

// -----------------------------------------------------------------------------
#define HANDLE_VERSION_OPTION()                                                \
  PARSE_STANDARD_OPTION()

// -----------------------------------------------------------------------------
#define HANDLE_UNKNOWN_OPTION()                                                \
  do {                                                                         \
    PrintHelp(EXECNAME);                                                       \
    std::cout << std::endl;                                                    \
    std::cout.flush();                                                         \
    std::cerr << "Error: Unknown option " << argv[OPTIDX] << std::endl;        \
    exit(1);                                                                   \
  } while (false)

// -----------------------------------------------------------------------------
#define HANDLE_STANDARD_OPTION()                                               \
  do {                                                                         \
    if      (HELP_OPTION     ) HANDLE_HELP_OPTION();                           \
    else if (VERSION_OPTION  ) HANDLE_VERSION_OPTION();                        \
    else if (STANDARD_OPTION ) PARSE_STANDARD_OPTION();                        \
  } while (false)

// -----------------------------------------------------------------------------
#define HANDLE_STANDARD_OR_UNKNOWN_OPTION()                                    \
  do {                                                                         \
    if      (HELP_OPTION     ) HANDLE_HELP_OPTION();                           \
    else if (VERSION_OPTION  ) HANDLE_VERSION_OPTION();                        \
    else if (STANDARD_OPTION ) PARSE_STANDARD_OPTION();                        \
    else                       HANDLE_UNKNOWN_OPTION();                        \
  } while (false)

// -----------------------------------------------------------------------------
#define HANDLE_COMMON_OPTION()                                                 \
  do {                                                                         \
    if      (HELP_OPTION     ) HANDLE_HELP_OPTION();                           \
    else if (VERSION_OPTION  ) HANDLE_VERSION_OPTION();                        \
    else if (STANDARD_OPTION ) PARSE_STANDARD_OPTION();                        \
    else if (PARALLEL_OPTION ) PARSE_PARALLEL_OPTION();                        \
    else if (PROFILING_OPTION) PARSE_PROFILING_OPTION();                       \
    else if (TERMINAL_OPTION ) PARSE_TERMINAL_OPTION();                        \
  } while (false)

// -----------------------------------------------------------------------------
#define HANDLE_COMMON_OR_UNKNOWN_OPTION()                                      \
  do {                                                                         \
    if      (HELP_OPTION     ) HANDLE_HELP_OPTION();                           \
    else if (VERSION_OPTION  ) HANDLE_VERSION_OPTION();                        \
    else if (STANDARD_OPTION ) PARSE_STANDARD_OPTION();                        \
    else if (PARALLEL_OPTION ) PARSE_PARALLEL_OPTION();                        \
    else if (PROFILING_OPTION) PARSE_PROFILING_OPTION();                       \
    else if (TERMINAL_OPTION ) PARSE_TERMINAL_OPTION();                        \
    else                       HANDLE_UNKNOWN_OPTION();                        \
  } while (false)

// -----------------------------------------------------------------------------
#define HANDLE_HELP_OR_VERSION()                                               \
  do {                                                                         \
    for (ALL_ARGUMENTS) {                                                      \
      if      (HELP_OPTION   ) HANDLE_HELP_OPTION();                           \
      else if (VERSION_OPTION) HANDLE_VERSION_OPTION();                        \
    }                                                                          \
  } while (false)

// =============================================================================
// Private global variables/functions used by command-line parsing macros
// =============================================================================

extern int         _posargc;                // Number of required positional arguments
extern int         _numposarg;              // Number of positional arguments
extern bool        _discard_parsed_posargs; // Whether to discard or keep parsed positional arguments
extern bool        _discard_parsed_options; // Whether to discard or keep parsed options
extern const char *_option;                 // Current option being parsed

int _GetNumberOfPositionalArguments(int, char *[]);

void  _DiscardArgument      (int &, int &, char *[]);
bool  _IsOption             (int &, int &, char *[], const char * = NULL);
bool  _IsArgument           (int,   int &, char *[]);
char *_GetPositionalArgument(int,   int &, char *[]);
char *_GetOptionArgument    (int &, int &, char *[]);


#endif
