#!/usr/bin/env python

from subprocess import check_call, check_output, STDOUT
from sys import argv, exit
from os.path import basename
from vtk import vtkPolyDataReader, vtkPolyDataWriter
import re


# ----------------------------------------------------------------------------------------------------------------------
def mris2vtk(input_name, output_name):
  """
  Converts a FreeSurfer surface to VTK file format
  """
  # convert surface to VTK format
  check_call(['mris_convert', input_name, output_name])
  # get surface RAS translation
  out = check_output(["mris_info", input_name], stderr=STDOUT)
  m = re.search("c_\(ras\)\s:\s\((-?\d+\.\d+),\s(-?\d+\.\d+),\s(-?\d+\.\d+)\)", out)
  if m is None: raise RuntimeError('Could not find c_(ras) coordinates in mris_info output!')
  tx = float(m.group(1))
  ty = float(m.group(2))
  tz = float(m.group(3))
  # transform vertex positions to scanner RAS of orig.mgz
  reader = vtkPolyDataReader()
  reader.SetFileName(output_name)
  reader.Update()
  surface = reader.GetOutput()
  points = surface.GetPoints()
  for i in range(points.GetNumberOfPoints()):
    x, y, z = points.GetPoint(i)
    points.SetPoint(i, x + tx, y + ty, z + tz)
  surface.SetPoints(points)
  writer = vtkPolyDataWriter()
  writer.SetFileName(output_name)
  writer.SetInput(surface)
  writer.Write()

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
  if len(argv) != 3:
    print "usage: " + basename(argv[0]) + " <surf> <output>"
    exit(1)
  mris2vtk(argv[1], argv[2])
