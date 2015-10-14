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

#include <irtkCxxLib.h>


// ========================================================================
// Swap bytes
// ========================================================================

// ------------------------------------------------------------------------
void swap16(char *a, char *b, long n)
{
  int i;
  char c;

  for (i = 0; i < n * 2; i += 2) {
    c = a[i];
    a[i] = b[i+1];
    b[i+1] = c;
  }
}

// ------------------------------------------------------------------------
void swap32(char *a, char *b, long n)
{
  int i;
  char c;

  for (i = 0; i < n * 4; i += 4) {
    c = a[i];
    a[i] = b[i+3];
    b[i+3] = c;
    c = a[i+1];
    a[i+1] = b[i+2];
    b[i+2] = c;
  }
}

// ------------------------------------------------------------------------
void swap64(char *a, char *b, long n)
{
  int i;
  char c;

  for (i = 0; i < n * 8; i += 8) {
    c = a[i];
    a[i] = b[i+7];
    b[i+7] = c;
    c = a[i+1];
    a[i+1] = b[i+6];
    b[i+6] = c;
    c = a[i+2];
    a[i+2] = b[i+5];
    b[i+5] = c;
    c = a[i+3];
    a[i+3] = b[i+4];
    b[i+4] = c;
  }
}

// ========================================================================
// Read from file
// ========================================================================

#define MAX_LINE 255

// ------------------------------------------------------------------------
int ReadInt(ifstream &in)
{
  char c;
  char *string, *s;
  int data;

  string = new char[MAX_LINE];
  s = string;
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadInt: No valid line format\n";
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  data = atoi(string);
  delete s;
  return (data);
}

// ------------------------------------------------------------------------
float ReadFloat(ifstream &in)
{
  char c;
  char *string, *s;
  float data;

  string = new char[MAX_LINE];
  s = string;
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadFloat: No valid line format\n";
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  data = atof(string);
  delete s;
  return (data);
}

// ------------------------------------------------------------------------
char *ReadString(ifstream &in)
{
  char c;
  char *string;

  string = new char[MAX_LINE];
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadString: No valid line format\n";
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  return (string);
}

// ------------------------------------------------------------------------
vector<string> Split(string s, const char *d, int n)
{
  const size_t delimiter_length = strlen(d);
  vector<string> parts;
  if (n <= 0) {
    n = abs(n);
    size_t pos;
    while (n == 0 || parts.size() < static_cast<size_t>(n)) {
      pos = s.rfind(d);
      if (pos == string::npos) {
        parts.push_back(s);
        break;
      }
      parts.push_back(s.substr(pos + delimiter_length));
      s.erase(pos);
    }
    reverse(parts.begin(), parts.end());
  } else {
    size_t start = 0;
    size_t end   = string::npos;
    while (start < s.length() && (n == 0 || parts.size() < static_cast<size_t>(n))) {
      end = s.find(d, start);
      parts.push_back(s.substr(start, end));
      if (end == string::npos) break;
      start = end + delimiter_length;
    }
  }
  return parts;
}
