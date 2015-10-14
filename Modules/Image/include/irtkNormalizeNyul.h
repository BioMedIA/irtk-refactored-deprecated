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
#ifndef IRTKNORMALIZENYUL_H_
#define IRTKNORMALIZENYUL_H_

#include <irtkImage.h>
//#include <irtkRegistration.h>

class irtkNormalizeNyul{

private:
	irtkRealImage _target;
	irtkRealImage _source;
	int _source_padding;
	int _target_padding;



public:
	irtkNormalizeNyul(irtkRealImage source, irtkRealImage target);
	void SetMask(irtkRealImage source_mask, irtkRealImage target_mask);
	void SetPadding(int source_padding, int target_padding);
	void Run();
	irtkRealImage GetOutput();
	irtkRealImage GetTarget();
};

#endif
