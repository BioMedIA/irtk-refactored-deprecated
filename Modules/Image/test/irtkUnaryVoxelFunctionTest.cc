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

#include <gtest/gtest.h>

#include <irtkImage.h>

#include <irtkVoxelFunction.h>
using namespace irtkImageDomain;

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
TEST(irtkUnaryVoxelFunction, GetMin)
{
  irtkGreyImage image;
  // Empty image
  {
    irtkUnaryVoxelFunction::GetMin min;
    ForEachVoxel(image, min);
    EXPECT_TRUE(IsNaN(min.GetMinAsDouble()));
  }
  // 2D image
  image.Initialize(11, 7);
  image.Put(2, 6,   7);
  image.Put(5, 3, -42);
  image.Put(2, 3, -41);
  {
    irtkUnaryVoxelFunction::GetMin min;
    ForEachVoxel(image, min);
    EXPECT_EQ(-42.0, min.GetMinAsDouble());
  }
  {
    irtkUnaryVoxelFunction::GetMin min;
    ParallelForEachVoxel(image, min);
    EXPECT_EQ(-42.0, min.GetMinAsDouble());
  }
  {
    irtkUnaryVoxelFunction::GetMin min;
    ForEachVoxelIf<Foreground>(image, min);
    EXPECT_EQ(-42.0, min.GetMinAsDouble());
  }
  {
    irtkUnaryVoxelFunction::GetMin min;
    ParallelForEachVoxelIf<Foreground>(image, min);
    EXPECT_EQ(-42.0, min.GetMinAsDouble());
  }
  image.PutBackgroundValueAsDouble(.0);
  {
    irtkUnaryVoxelFunction::GetMin min;
    ForEachVoxelIf<Foreground>(image, min);
    EXPECT_EQ(-42.0, min.GetMinAsDouble());
    min.Reset();
    ForEachVoxelIf<AboveBackgroundLevel>(image, min);
    EXPECT_EQ(7.0, min.GetMinAsDouble());
  }
  {
    irtkUnaryVoxelFunction::GetMin min;
    ParallelForEachVoxelIf<Foreground>(image, min);
    EXPECT_EQ(-42.0, min.GetMinAsDouble());
    min.Reset();
    ParallelForEachVoxelIf<AboveBackgroundLevel>(image, min);
    EXPECT_EQ(7.0, min.GetMinAsDouble());
  }
  image.ClearBackgroundValue();
  // 3D image
  /*
  image.Initialize(11, 7, 5);
  image.Put(2, 3, 0,   7);
  image.Put(5, 3, 3, -42);
  image.Put(2, 3, 2, -41);
  // 4D image
  image.Initialize(11, 7, 5, 3);
  */
}

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
