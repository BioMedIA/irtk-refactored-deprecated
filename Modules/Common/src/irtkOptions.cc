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
#include <irtkVersionInfo.h>


// =============================================================================
// Global standard options
// =============================================================================

const irtkVersion current_version(IRTK_VERSION_MAJOR,
                                  IRTK_VERSION_MINOR,
                                  IRTK_VERSION_PATCH);

// Default: "Emulate" latest version
irtkVersion version = current_version;

// Default: No verbose output messages
int verbose = 0;

// Default: No debug output
int debug = 0;

// =============================================================================
// Help
// =============================================================================

int         _posargc                = 0;
int         _numposarg              = -1;
bool        _discard_parsed_posargs = false;
bool        _discard_parsed_options = false;
const char *_option                 = NULL;

// -----------------------------------------------------------------------------
void PrintRevision(ostream &out)
{
  if (IRTK_REVISION) out << IRTK_REVISION;
  else               out << current_version << endl;
}

// -----------------------------------------------------------------------------
void PrintVersion(ostream &out, const char *name)
{
  if (name) out << name << " ";
  if (current_version) {
    out << current_version;
    if (version < current_version) out << ", emulating version " << version;
  }
  out << " (";
  if (IRTK_REVISION) out << "rev " << IRTK_REVISION << ", ";
  out << "built on " << (__DATE__);
  out << ")";
  out << endl;
}

// -----------------------------------------------------------------------------
bool IsStandardOption(const char *arg)
{
  _option = NULL;
  if      (strcmp(arg, "-v")        == 0) _option = "-v";
  else if (strcmp(arg, "-verbose")  == 0) _option = "-verbose";
  else if (strcmp(arg, "-debug")    == 0) _option = "-debug";
  else if (strcmp(arg, "-revision") == 0) _option = "-revision";
  else if (strcmp(arg, "-version")  == 0 || strcmp(arg, "--version") == 0) _option = "-version";
  return (_option != NULL);
}

// -----------------------------------------------------------------------------
void ParseStandardOption(int &OPTIDX, int &argc, char *argv[])
{
  if (OPTION("-v") || OPTION("-verbose")) {
    if (HAS_ARGUMENT) verbose  = atoi(ARGUMENT);
    else              verbose += 1;
  } else if (OPTION("-debug")) {
    if (HAS_ARGUMENT) debug  = atoi(ARGUMENT);
    else              debug += 1;
  } else if (OPTION("-version") || OPTION("--version")) {
    if (HAS_ARGUMENT) {
      if (!FromString(ARGUMENT, version) || version > current_version) {
        cerr << "Invalid [-]-version argument" << endl;
        exit(1);
      }
    } else {
      PrintVersion(cout, EXECNAME);
      exit(0);
    }
  } else if (OPTION("-revision")) {
    PrintRevision(cout);
    exit(0);
  }
}

// -----------------------------------------------------------------------------
void PrintStandardOptions(ostream &out)
{
  out << endl;
  out << "Standard options:" << endl;
  out << "  -v -verbose [n]              Increase/Set verbosity of output messages. (default: " << verbose << ")" << endl;
  out << "  -debug [level]               Increase/Set debug level for output of intermediate results. (default: " << debug << ")" << endl;
  out << "  -[-]version [major.minor]    Print version and exit or set version to emulate." << endl;
  out << "  -revision                    Print revision (or version) number only and exit." << endl;
  out << "  -h -[-]help                  Print help and exit." << endl;
}

// -----------------------------------------------------------------------------
void PrintCommonOptions(ostream &out)
{
  PrintStandardOptions (out);
  PrintTerminalOptions (out);
  PrintParallelOptions (out);
  PrintProfilingOptions(out);
}

// =============================================================================
// Private functions used by command-line parsing macros
// =============================================================================

// -----------------------------------------------------------------------------
int _GetNumberOfPositionalArguments(int argc, char *argv[])
{
  int n = 1;
  while (n < argc && argv[n][0] != '-') n++;
  return n - 1;
}

// -----------------------------------------------------------------------------
void _DiscardArgument(int &i, int &argc, char *argv[])
{
  for (int j = i; j < argc; j++) argv[j] = argv[j+1];
  argv[argc--] = NULL;
  i--;
}

// -----------------------------------------------------------------------------
bool _IsOption(int &i, int &argc, char *argv[], const char *opt)
{
  if ((opt == NULL && argv[i][0] == '-') || (opt != NULL && strcmp(argv[i], opt) == 0)) {
    if (_discard_parsed_options) _DiscardArgument(i, argc, argv);
    _option = argv[i];
    return true;
  } else {
    _option = NULL;
    return false;
  }
}

// -----------------------------------------------------------------------------
bool _IsArgument(int i, int &argc, char *argv[])
{
  return (i+1 < argc && argv[i+1][0] != '-');
}

// -----------------------------------------------------------------------------
char *_GetPositionalArgument(int i, int &argc, char *argv[])
{
  if (i >= argc) {
    cerr << "Error: Not all required arguments specified!" << endl;
    exit(1);
  }
  if (_posargc < i) _posargc = i;
  char *arg = argv[i];
  return arg;
}

// -----------------------------------------------------------------------------
char *_GetOptionArgument(int &i, int &argc, char *argv[])
{
  if (_option) {
    i++;
    if (i >= argc) {
      cerr << "Error: Option";
      if (_option) cerr << " " << _option;
      cerr << " requires more arguments than specified!" << endl;
      exit(1);
    }
  } else {
    if (i >= argc) {
      cerr << "Error: Not all required arguments specified!" << endl;
      exit(1);
    }
    if (_posargc < i) _posargc = i;
  }
  char *arg = argv[i];
  if (_discard_parsed_options) _DiscardArgument(i, argc, argv);
  return arg;
}
