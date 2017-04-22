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
 * Neither the name of the copyright holder nor the names of its
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
#include <CameraModels/PinholeDistance.hpp>

/**
 * Pinhole, using distance for inverse.
 */
#include <CameraModels/PinholeDistanceDistorted.hpp> 

/**
 * Pinhole no distortions, inverse in depth/planes (not distance).
 */
#include <CameraModels/PinholeDisparity.hpp> 

/**
 * Pinhole, inverse in image planes (not distance).
 */
#include <CameraModels/PinholeDisparityDistorted.hpp> 

/**
 * Gayer-Baretto-Mei Model, no distortions.
 */
#include <CameraModels/Generic.hpp>

/**
 * Gayer-Baretto-Mei Model.
 */
#include <CameraModels/GenericDistorted.hpp> 

/**
 * Simple Spherical Panorama camera model.
 */
#include <CameraModels/Spherical.hpp>

/**
 * Spherical Panorama camera model as in PovRay.
 */
#include <CameraModels/SphericalPovRay.hpp>

/**
 * OpenCV 3.0 Fisheye model.
 */
#include <CameraModels/FisheyeDistorted.hpp> 

/**
 * OpenCV 3.0 Fisheye model, no distortions.
 */
#include <CameraModels/Fisheye.hpp> 

/**
 * Pinhole, inverse in image planes (not distance), modified Brown-Conrady distortions (cf Intel librealsense).
 * @note not invertible.
 */
#include <CameraModels/PinholeDisparityBrownConrady.hpp> 

namespace cammod
{

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDistance>
{
    template<typename T> using ModelT = PinholeDistance<T>;
    static constexpr auto Name = "PinholeDistance";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDistanceDistorted>
{
    template<typename T> using ModelT = PinholeDistanceDistorted<T>;
    static constexpr auto Name = "PinholeDistanceDistorted";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDisparity>
{
    template<typename T> using ModelT = PinholeDisparity<T>;
    static constexpr auto Name = "PinholeDisparity";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDisparityDistorted>
{
    template<typename T> using ModelT = PinholeDisparityDistorted<T>;
    static constexpr auto Name = "PinholeDisparityDistorted";
};

template<> struct CameraModelToTypeAndName<CameraModelType::Generic>
{
    template<typename T> using ModelT = Generic<T>;
    static constexpr auto Name = "Generic";
};

template<> struct CameraModelToTypeAndName<CameraModelType::GenericDistorted>
{
    template<typename T> using ModelT = GenericDistorted<T>;
    static constexpr auto Name = "GenericDistorted";
};

template<> struct CameraModelToTypeAndName<CameraModelType::Spherical>
{
    template<typename T> using ModelT = Spherical<T>;
    static constexpr auto Name = "Spherical";
};

template<> struct CameraModelToTypeAndName<CameraModelType::SphericalPovRay>
{
    template<typename T> using ModelT = SphericalPovRay<T>;
    static constexpr auto Name = "SphericalPovRay";
};

template<> struct CameraModelToTypeAndName<CameraModelType::Fisheye>
{
    template<typename T> using ModelT = Fisheye<T>;
    static constexpr auto Name = "Fisheye";
};

template<> struct CameraModelToTypeAndName<CameraModelType::FisheyeDistorted>
{
    template<typename T> using ModelT = FisheyeDistorted<T>;
    static constexpr auto Name = "FisheyeDistorted";
};

template<> struct CameraModelToTypeAndName<CameraModelType::PinholeDisparityBrownConrady>
{
    template<typename T> using ModelT = PinholeDisparityBrownConrady<T>;
    static constexpr auto Name = "PinholeDisparityBrownConrady";
};

}

#endif // CAMERA_MODELS_HPP
