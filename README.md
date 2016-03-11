CameraModels
=============

This is a small header only library of templated C++11 classes implementing 
various camera models.

Based on Eigen and Sophus for algebra and 3D transformations.

## Features
Each class provides both normal and static methods for the following functions:

* `forward (3d_point) -> pixel`,
* `forward (pose, 3d_point) -> pixel`,
* `inverse (pixel) -> 3d_point`, either unit vector or on image plane Z=1,
* `inverseAtDistance (pixel, distance) -> 3d_point`,
* `inverseAtDistance (pose, pixel, distance) -> 3d_point`,
* `twoFrameProject (pose1, pixel1, pose2) -> 3d_point`,
* `worldToCamera(pose, 3d_point) -> 3d_point`,
* `cameraToWorld(pose, 3d_point) -> 3d_point`,
* `pixelValid(pixel) -> bool`,
* `cast<NewScalarType>()`,
* `resizeViewport(dimensions)`.

##### Type safety
Each model is templated with the scalar type (e.g. float, double). Statically
typed models can work across CUDA host-device boundary provided the underlying 
infrastructure is CUDA ready (Eigen + Sophus). Additionally, each model can be 
wrapped into _Eigen::Map_ (useful when used inside Ceres-Solver cost functions).

While each model is typed with the default scalar type, most of the methods can be
retyped with a different type. For example, when doing camera calibration with 
Ceres-Solver automatic differentiation typing the entire model with _ceres::Jet_
is necessary, while when solving Bundle-Adjustment problem the camera model class
can stay typed with double, while only the _forward<ceres::Jet>(...)_ call 
being typed with _ceres::Jet_ to carry out the derivation.

##### Conventions
Poses are _Sophus::SE3Group_ Lie algebra transformations, points and pixels are Eigen
3D and 2D vectors (see [CameraModelHelpers.hpp](CameraModelHelpers.hpp) for conventions).

##### Calibration
Some models do not support calibration, as it wasn't tested. Each model provides
the following constants:
* NumParameters - total number of parameters, including miscellanea like width & height,
* ParametersToOptimize - first N parameters, only the ones to optimize, 
* CalibrationSupported - true/false, if the model supports calibration.

Most models were tested with Ceres-Solver based calibration employing automatic differentiation.

#####  Dynamic Polymorphism
If run time polymorphism is needed for whatever reason one might use _CameraInterface_ class
and wrap statically typed model with _CameraFromCRTP_ that exposes virtual functions.

##### Misc. features
Each camera model has _std::ostream_ operator. Optionally, Boost Serialization or
Cereal serialization is supported (see [CameraModelHelpers.hpp](CameraModelHelpers.hpp) for more info)
Camera models with distortions implement _getIdeal()_ method that returns the
corresponding camera model class with distortions omitted. Note: this is a very crude
way of obtaining the ideal camera model to work on undistorted images.

## Models supported

* Pinhole - classical pinhole camera model, however inverts with distance,
* PinholeDistorted - as above, but with radial-tangential distortions,
* IdealGeneric - Gayer/Baretto/Mei generic central camera model,
* FullGeneric - as above, but with radial-tangential distortions,
* Spherical - simple spherical panorama,
* SphericalPovRay - as above, but following PovRay conventions,
* IdealFisheye - fisheye camera model as seen in OpenCV 3.0,
* Fisheye - as above, including extra distortion coefficients,
* PinholeDisparity - classical pinhole camera model, inverse on the image planes,
* PinholeDisparityDistorted - as above, but with radial-tangential distortions,
* PinholeDisparityBrownConrady - special case, non invertible. See Intel's librealsense.

## TODO
* Much more testing, not just simple forward/inverse checks,
* Additional projection-unprojection checking functions, tests for chirality etc.,
* Implement camera rig type to represent multi-camera systems.

## License

_Copyright © `2015`, `Robert Lukierski`_  
_All rights reserved._

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.