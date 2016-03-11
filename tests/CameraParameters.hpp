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
 * Sensible values to initialize camera models to test.
 * ****************************************************************************
 */

#ifndef CAMERA_PARAMETERS_HPP
#define CAMERA_PARAMETERS_HPP

#include <CameraModels.hpp>


template<typename ModelT>
struct CameraParameters
{
};

template<> struct CameraParameters<camera::PinholeCameraModel<float>>
{
    typedef camera::PinholeCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeCameraModel<double>>
{
    typedef camera::PinholeCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDisparityCameraModel<float>>
{
    typedef camera::PinholeDisparityCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDisparityCameraModel<double>>
{
    typedef camera::PinholeDisparityCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDistortedCameraModel<float>>
{
    typedef camera::PinholeDistortedCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, -2.460504e-01f, 9.947797e-02f, -1.070869e-03f, -5.855505e-06f, 0.000000e+00f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDistortedCameraModel<double>>
{
    typedef camera::PinholeDistortedCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, -2.460504e-01, 9.947797e-02, -1.070869e-03, -5.855505e-06, 0.000000e+00, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDisparityDistortedCameraModel<float>>
{
    typedef camera::PinholeDisparityDistortedCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, -2.460504e-01f, 9.947797e-02f, -1.070869e-03f, -5.855505e-06f, 0.000000e+00f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDisparityDistortedCameraModel<double>>
{
    typedef camera::PinholeDisparityDistortedCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, -2.460504e-01, 9.947797e-02, -1.070869e-03, -5.855505e-06, 0.000000e+00, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDisparityBrownConradyCameraModel<float>>
{
    typedef camera::PinholeDisparityBrownConradyCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 640;
    static constexpr std::size_t DefaultHeight = 512;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(381.752f, 381.866f, 323.454f, 268.756f, -0.358691f, 0.161264f, -0.037037f, 0.000249145f, 0.000462129f, 0.0f, (float)DefaultWidth, (float)DefaultHeight);
    }
};

template<> struct CameraParameters<camera::PinholeDisparityBrownConradyCameraModel<double>>
{
    typedef camera::PinholeDisparityBrownConradyCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 640;
    static constexpr std::size_t DefaultHeight = 512;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(381.752, 381.866, 323.454, 268.756, -0.358691, 0.161264, -0.037037, 0.000249145, 0.000462129, 0.0, (double)DefaultWidth, (double)DefaultHeight);
    }
};

template<> struct CameraParameters<camera::FullGenericCameraModel<float>>
{
    typedef camera::FullGenericCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, 0.75f, -2.460504e-01f, 9.947797e-02f, -1.070869e-03f, -5.855505e-06f, 0.000000e+00f, DefaultWidth, DefaultHeight, 180, 350);
    }
};

template<> struct CameraParameters<camera::FullGenericCameraModel<double>>
{
    typedef camera::FullGenericCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, 0.75, -2.460504e-01, 9.947797e-02, -1.070869e-03, -5.855505e-06, 0.000000e+00, DefaultWidth, DefaultHeight, 180, 350);
    }
};

template<> struct CameraParameters<camera::IdealGenericCameraModel<float>>
{
    typedef camera::IdealGenericCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, 0.75f, DefaultWidth, DefaultHeight, 180, 350);
    }
};

template<> struct CameraParameters<camera::IdealGenericCameraModel<double>>
{
    typedef camera::IdealGenericCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, 0.75, DefaultWidth, DefaultHeight, 180, 350);
    }
};

static constexpr float FisheyeRadius = 500.0f;
static constexpr float IdealFisheyeRadius = 495.0f;

template<> struct CameraParameters<camera::FisheyeCameraModel<float>>
{
    typedef camera::FisheyeCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866f, 315.678f, 764.713f, 580.514f, 0.019788f, -0.00758584f, 0.00235877f, -0.000606299f, -0.000586391f, DefaultWidth, DefaultHeight, FisheyeRadius);
    }
};

template<> struct CameraParameters<camera::FisheyeCameraModel<double>>
{
    typedef camera::FisheyeCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866, 315.678, 764.713, 580.514, 0.019788, -0.00758584, 0.00235877, -0.000606299, -0.000586391, DefaultWidth, DefaultHeight, FisheyeRadius);
    }
};

template<> struct CameraParameters<camera::IdealFisheyeCameraModel<float>>
{
    typedef camera::IdealFisheyeCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866f, 315.678f, 764.713f, 580.514f, DefaultWidth, DefaultHeight, IdealFisheyeRadius);
    }
};

template<> struct CameraParameters<camera::IdealFisheyeCameraModel<double>>
{
    typedef camera::IdealFisheyeCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866, 315.678, 764.713, 580.514, DefaultWidth, DefaultHeight, IdealFisheyeRadius);
    }
};

template<> struct CameraParameters<camera::SphericalCameraModel<float>>
{
    typedef camera::SphericalCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0f, 256.0f, 0.895866f, 1.80258f);
    }
};

template<> struct CameraParameters<camera::SphericalCameraModel<double>>
{
    typedef camera::SphericalCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0, 256.0, 0.895866, 1.80258);
    }
};

template<> struct CameraParameters<camera::SphericalPovRayCameraModel<float>>
{
    typedef camera::SphericalPovRayCameraModel<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0f, 256.0f, 0.895866f, 1.80258f);
    }
};

template<> struct CameraParameters<camera::SphericalPovRayCameraModel<double>>
{
    typedef camera::SphericalPovRayCameraModel<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0, 256.0, 0.895866, 1.80258);
    }
};

#endif // CAMERA_PARAMETERS_HPP