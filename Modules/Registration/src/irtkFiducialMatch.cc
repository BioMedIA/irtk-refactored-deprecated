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

#include <irtkFiducialMatch.h>

#include <irtkPoint.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkFiducialMatch::irtkFiducialMatch()
{
}

// -----------------------------------------------------------------------------
irtkFiducialMatch
::irtkFiducialMatch(const irtkFiducialMatch &other)
:
  irtkPointCorrespondence(other)
{
}

// -----------------------------------------------------------------------------
irtkPointCorrespondence *irtkFiducialMatch::NewInstance() const
{
  return new irtkFiducialMatch(*this);
}

// -----------------------------------------------------------------------------
irtkFiducialMatch::~irtkFiducialMatch()
{
}

// -----------------------------------------------------------------------------
irtkFiducialMatch::TypeId irtkFiducialMatch::Type() const
{
  return FiducialMatch;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkFiducialMatch::Set(const char *name, const char *value)
{
  if (strcmp(name, "Correspondence map") == 0) {
    _CorrespondenceMap = value;
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
irtkParameterList irtkFiducialMatch::Parameter() const
{
  irtkParameterList params;
  Insert(params, "Correspondence map", _CorrespondenceMap);
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void irtkFiducialMatch::ValidateCorrespondenceMap(const irtkRegisteredPointSet *target,
                                                  const irtkRegisteredPointSet *source,
                                                  const vector<int>            &map,
                                                  const char                   *map_name) const
{
  const int m = static_cast<int>(target->NumberOfPoints());
  const int n = static_cast<int>(source->NumberOfPoints());
  if (map.size() < static_cast<size_t>(m)) {
    cerr << NameOfType() << "::Initialize: Correspondence map has no entries for target indices t >= " << map.size() << endl;
    if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
    exit(1);
  }
  if (map.size() > static_cast<size_t>(m)) {
    cerr << NameOfType() << "::Initialize: Correspondence map has additional entries for unused indices t >= " << m << endl;
    if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
  }
  for (int t = 0; t < m; ++t) {
    if (map[t] == -1) {
      cerr << NameOfType() << "::Initialize: Missing correspondence map entry for t=" << t << endl;
      if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
      exit(1);
    } else if (map[t] < 0 || map[t] >= n) {
      cerr << NameOfType() << "::Initialize: Invalid correspondence map entry: t=" << t << ", s=" << map[t] << endl;
      if (map_name && map_name[0] != '\0') cerr << "  Correspondence map input file: " << map_name << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void irtkFiducialMatch::InvertCorrespondenceMap(const irtkRegisteredPointSet *target,
                                                const irtkRegisteredPointSet *source,
                                                const vector<int>            &map,
                                                vector<int>                  &inv) const
{
  const int m = static_cast<int>(target->NumberOfPoints());
  const int n = static_cast<int>(source->NumberOfPoints());
  inv.clear();
  inv.resize(n, -1);
  for (int t = 0; t < m; ++t) {
    inv[map[t]] = t;
  }
}

// -----------------------------------------------------------------------------
void irtkFiducialMatch::Initialize()
{
  // Initialize base class
  irtkPointCorrespondence::Initialize();

  const int m = static_cast<int>(_Target->NumberOfPoints());
  const int n = static_cast<int>(_Source->NumberOfPoints());

  // Either read corresponding indices from input file
  if (!_CorrespondenceMap.empty()) {
    ifstream ifs(_CorrespondenceMap);
    if (!ifs.is_open()) {
      cerr << NameOfType()
           << "::Initialize: Failed to open correspondence input file: "
           << _CorrespondenceMap << endl;
      exit(1);
    }
    string line;
    int    t, s;
    _SourceIndex.resize(m, -1);
    while (getline(ifs, line)) {
      istringstream is(line);
      if (!(is >> t >> s)) {
        cerr << NameOfType()
             << "::Initialize: Failed to read correspondence map from file: "
             << _CorrespondenceMap << endl;
        exit(1);
      }
      if (t < 0 || t >= m) {
        cerr << NameOfType()
             << "::Initialize: Invalid target index in correspondence map: "
             << _CorrespondenceMap << endl;
        exit(1);
      }
      _SourceIndex[t] = s;
    }
    ifs.close();
  // Or ensure that both data sets contain same number of fiducial points
  } else {
    if (m != n) {
      cerr << NameOfType() << "::Initialize: Data sets must have same number of fiducial points" << endl;
      exit(1);
    }
    _SourceIndex.clear();
    _TargetIndex.clear();
  }

  if (!_SourceIndex.empty()) {
    // Check correspondence map
    ValidateCorrespondenceMap(_Target, _Source, _SourceIndex, _CorrespondenceMap.c_str());
    // Compute inverse correspondence map
    InvertCorrespondenceMap(_Target, _Source, _SourceIndex, _TargetIndex);
    // Check inverse correspondence map
    ValidateCorrespondenceMap(_Source, _Target, _TargetIndex);
  }
}

// -----------------------------------------------------------------------------
bool irtkFiducialMatch::GetInputTargetPoint(int i, irtkPoint &p) const
{
  if (_SourceSample) i = (*_SourceSample)[i];
  if (!_TargetIndex.empty()) i = _TargetIndex[i];
  _Target->GetInputPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkFiducialMatch::GetInputSourcePoint(int i, irtkPoint &p) const
{
  if (_TargetSample) i = (*_TargetSample)[i];
  if (!_SourceIndex.empty()) i = _SourceIndex[i];
  _Source->GetInputPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkFiducialMatch::GetTargetPoint(int i, irtkPoint &p) const
{
  if (_SourceSample) i = (*_SourceSample)[i];
  if (!_TargetIndex.empty()) i = _TargetIndex[i];
  _Target->GetPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkFiducialMatch::GetSourcePoint(int i, irtkPoint &p) const
{
  if (_TargetSample) i = (*_TargetSample)[i];
  if (!_SourceIndex.empty()) i = _SourceIndex[i];
  _Source->GetPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
int irtkFiducialMatch::GetTargetIndex(int i) const
{
  if (_SourceSample) i = (*_SourceSample)[i];
  if (!_TargetIndex.empty()) i = _TargetIndex[i];
  return i;
}

// -----------------------------------------------------------------------------
int irtkFiducialMatch::GetSourceIndex(int i) const
{
  if (_TargetSample) i = (*_TargetSample)[i];
  if (!_SourceIndex.empty()) i = _SourceIndex[i];
  return i;
}
