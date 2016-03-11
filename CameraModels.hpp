/**
 * ****************************************************************************
 * Copyright (c) 2015, Robert Lukierski.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * Neither the name of [project] nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * ****************************************************************************
 * Top Level Header.
 * ****************************************************************************
 */

#ifndef CAMERA_MODELS_HPP
#define CAMERA_MODELS_HPP

/**
 * Pinhole no distortions, using distance for inverse.
 */
#include <PinholeCameraModel.hpp>

/**
 * Pinhole, using distance for inverse.
 */
#include <PinholeDistortedCameraModel.hpp> 

/**
 * Gayer-Baretto-Mei Model, no distortions.
 */
#include <IdealGenericCameraModel.hpp>

/**
 * Gayer-Baretto-Mei Model.
 */
#include <FullGenericCameraModel.hpp> 

/**
 * Simple Spherical Panorama camera model.
 */
#include <SphericalCameraModel.hpp>

/**
 * Spherical Panorama camera model as in PovRay.
 */
#include <SphericalPovRayCameraModel.hpp>

/**
 * OpenCV 3.0 Fisheye model.
 */
#include <FisheyeCameraModel.hpp> 

/**
 * OpenCV 3.0 Fisheye model, no distortions.
 */
#include <IdealFisheyeCameraModel.hpp> 

/**
 * Pinhole no distortions, inverse in depth/planes (not distance).
 */
#include <PinholeDisparityCameraModel.hpp> 

/**
 * Pinhole, inverse in image planes (not distance).
 */
#include <PinholeDisparityDistortedCameraModel.hpp> 

/**
 * Pinhole, inverse in image planes (not distance), modified Brown-Conrady distortions (cf Intel librealsense).
 * @note not invertible.
 */
#include <PinholeDisparityBrownConrady.hpp> 

namespace camera
{

template<> struct CameraModelToTypeAndName<CameraModelType::Pinhole>
{
    template<typename T> using ModelT = PinholeCameraModel<T>;
    static constexpr auto Name = "Pinhole";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDistorted>
{
    template<typename T> using ModelT = PinholeDistortedCameraModel<T>;
    static constexpr auto Name = "PinholeDistorted";
};

template<> struct CameraModelToTypeAndName<CameraModelType::IdealGeneric>
{
    template<typename T> using ModelT = IdealGenericCameraModel<T>;
    static constexpr auto Name = "IdealGeneric";
};

template<> struct CameraModelToTypeAndName<CameraModelType::FullGeneric>
{
    template<typename T> using ModelT = FullGenericCameraModel<T>;
    static constexpr auto Name = "FullGeneric";
};

template<> struct CameraModelToTypeAndName<CameraModelType::Spherical>
{
    template<typename T> using ModelT = SphericalCameraModel<T>;
    static constexpr auto Name = "Spherical";
};

template<> struct CameraModelToTypeAndName<CameraModelType::SphericalPovRay>
{
    template<typename T> using ModelT = SphericalPovRayCameraModel<T>;
    static constexpr auto Name = "SphericalPovRay";
};

template<> struct CameraModelToTypeAndName<CameraModelType::Fisheye>
{
    template<typename T> using ModelT = FisheyeCameraModel<T>;
    static constexpr auto Name = "Fisheye";
};

template<> struct CameraModelToTypeAndName<CameraModelType::IdealFisheye>
{
    template<typename T> using ModelT = IdealFisheyeCameraModel<T>;
    static constexpr auto Name = "IdealFisheye";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDisparity>
{
    template<typename T> using ModelT = PinholeDisparityCameraModel<T>;
    static constexpr auto Name = "PinholeDisparity";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDisparityDistorted>
{
    template<typename T> using ModelT = PinholeDisparityDistortedCameraModel<T>;
    static constexpr auto Name = "PinholeDisparityDistorted";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDisparityBrownConrady>
{
    template<typename T> using ModelT = PinholeDisparityBrownConradyCameraModel<T>;
    static constexpr auto Name = "PinholeDisparityBrownConrady";
};

}

#endif // CAMERA_MODELS_HPP