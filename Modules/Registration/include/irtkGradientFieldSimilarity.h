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

#ifndef _IRTKGRADIENTFIELDSIMILARITY_H
#define _IRTKGRADIENTFIELDSIMILARITY_H

#include <irtkImageSimilarity.h>


/**
 * Base class for gradient field similarity measures
 *
 * Subclasses of this image similarity measure evaluate similarity of two images
 * based on their intensity gradient.
 */

class irtkGradientFieldSimilarity : public irtkImageSimilarity
{
  irtkObjectMacro(irtkGradientFieldSimilarity);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Noise parameter for target image
  irtkPublicAttributeMacro(bool, IgnoreJacobianGradientWrtDOFs);

  /// Transformed gradient of the target image
  /// Used only if the target image is being transformed
  irtkAttributeMacro(GradientImageType, TargetTransformedGradient);

  /// Transformed gradient of the source image
  /// Used only if the source image is being transformed
  irtkAttributeMacro(GradientImageType, SourceTransformedGradient);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  irtkGradientFieldSimilarity(const char * = "", double = 1.0);

  /// Copy constructor
  irtkGradientFieldSimilarity(const irtkGradientFieldSimilarity &);

  /// Destructor
  virtual ~irtkGradientFieldSimilarity();

  // ---------------------------------------------------------------------------
  // Parameters

public:

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Import other overloads
  using irtkImageSimilarity::Parameter;

  /// Get parameter name/value map
  virtual irtkParameterList Parameter() const;

protected:

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure once input and parameters have been set
  /// \param[in] domain Image domain on which the similarity is evaluated.
  virtual void InitializeInput(const irtkImageAttributes &domain);

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Convert non-parametric similarity gradient into gradient
  /// w.r.t transformation parameters
  ///
  /// This function calls irtkTransformation::ParametricGradient of the
  /// transformation to apply the chain rule in order to obtain the similarity
  /// gradient w.r.t the transformation parameters. It adds the weighted gradient
  /// to the final registration energy gradient.
  ///
  /// \param[in]    image       Transformed image.
  /// \param[in]    np_gradient Voxel-wise non-parametric gradient.
  /// \param[inout] gradient    Gradient to which the computed parametric gradient
  ///                           is added, after multiplication by the given \p weight.
  /// \param[in]    weight      Weight of image similarity.
  virtual void ParametricGradient(const irtkRegisteredImage *image,
                                  GradientImageType         *np_gradient,
                                  double                    *gradient,
                                  double                     weight);

  /// Update moving input image(s) and internal state of similarity measure
  virtual void Update(bool = true);

  /// Reorient transformed image gradient according to dI(y)/dy * dy/dx
  void ReorientGradient(irtkRegisteredImage *, bool = false);

  /// Multiply similarity gradient by 2nd order derivatives of transformed image
  ///
  /// \param[in]     image    Transformed image
  /// \param[in,out] gradient Input must be the gradient of the image similarity
  ///                         w.r.t. the transformed \p image gradient. Output is
  ///                         the voxel-wise gradient of the similarity w.r.t. T(x).
  void MultiplyByImageHessian(const irtkRegisteredImage *image,
                              GradientImageType         *gradient);

};


#endif
