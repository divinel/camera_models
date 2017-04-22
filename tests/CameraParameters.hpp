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
 * Sensible values to initialize camera models to test.
 * ****************************************************************************
 */

#ifndef CAMERA_PARAMETERS_HPP
#define CAMERA_PARAMETERS_HPP

#include <CameraModels/CameraModels.hpp>


template<typename ModelT>
struct CameraParameters
{
};

template<> struct CameraParameters<cammod::PinholeDistance<float>>
{
    typedef cammod::PinholeDistance<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDistance<double>>
{
    typedef cammod::PinholeDistance<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDisparity<float>>
{
    typedef cammod::PinholeDisparity<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDisparity<double>>
{
    typedef cammod::PinholeDisparity<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDistanceDistorted<float>>
{
    typedef cammod::PinholeDistanceDistorted<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, -2.460504e-01f, 9.947797e-02f, -1.070869e-03f, -5.855505e-06f, 0.000000e+00f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDistanceDistorted<double>>
{
    typedef cammod::PinholeDistanceDistorted<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, -2.460504e-01, 9.947797e-02, -1.070869e-03, -5.855505e-06, 0.000000e+00, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDisparityDistorted<float>>
{
    typedef cammod::PinholeDisparityDistorted<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, -2.460504e-01f, 9.947797e-02f, -1.070869e-03f, -5.855505e-06f, 0.000000e+00f, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDisparityDistorted<double>>
{
    typedef cammod::PinholeDisparityDistorted<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, -2.460504e-01, 9.947797e-02, -1.070869e-03, -5.855505e-06, 0.000000e+00, DefaultWidth, DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDisparityBrownConrady<float>>
{
    typedef cammod::PinholeDisparityBrownConrady<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 640;
    static constexpr std::size_t DefaultHeight = 512;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(381.752f, 381.866f, 323.454f, 268.756f, -0.358691f, 0.161264f, -0.037037f, 0.000249145f, 0.000462129f, 0.0f, (float)DefaultWidth, (float)DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::PinholeDisparityBrownConrady<double>>
{
    typedef cammod::PinholeDisparityBrownConrady<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 640;
    static constexpr std::size_t DefaultHeight = 512;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(381.752, 381.866, 323.454, 268.756, -0.358691, 0.161264, -0.037037, 0.000249145, 0.000462129, 0.0, (double)DefaultWidth, (double)DefaultHeight);
    }
};

template<> struct CameraParameters<cammod::GenericDistorted<float>>
{
    typedef cammod::GenericDistorted<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, 0.75f, -2.460504e-01f, 9.947797e-02f, -1.070869e-03f, -5.855505e-06f, 0.000000e+00f, DefaultWidth, DefaultHeight, 180, 350);
    }
};

template<> struct CameraParameters<cammod::GenericDistorted<double>>
{
    typedef cammod::GenericDistorted<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, 0.75, -2.460504e-01, 9.947797e-02, -1.070869e-03, -5.855505e-06, 0.000000e+00, DefaultWidth, DefaultHeight, 180, 350);
    }
};

template<> struct CameraParameters<cammod::Generic<float>>
{
    typedef cammod::Generic<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697f, 727.05397f, 649.83487f, 483.03012f, 0.75f, DefaultWidth, DefaultHeight, 180, 350);
    }
};

template<> struct CameraParameters<cammod::Generic<double>>
{
    typedef cammod::Generic<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1280;
    static constexpr std::size_t DefaultHeight = 960;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(746.96697, 727.05397, 649.83487, 483.03012, 0.75, DefaultWidth, DefaultHeight, 180, 350);
    }
};

static constexpr float FisheyeRadius = 500.0f;
static constexpr float IdealFisheyeRadius = 495.0f;

template<> struct CameraParameters<cammod::FisheyeDistorted<float>>
{
    typedef cammod::FisheyeDistorted<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866f, 315.678f, 764.713f, 580.514f, 0.019788f, -0.00758584f, 0.00235877f, -0.000606299f, -0.000586391f, DefaultWidth, DefaultHeight, FisheyeRadius);
    }
};

template<> struct CameraParameters<cammod::FisheyeDistorted<double>>
{
    typedef cammod::FisheyeDistorted<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866, 315.678, 764.713, 580.514, 0.019788, -0.00758584, 0.00235877, -0.000606299, -0.000586391, DefaultWidth, DefaultHeight, FisheyeRadius);
    }
};

template<> struct CameraParameters<cammod::Fisheye<float>>
{
    typedef cammod::Fisheye<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866f, 315.678f, 764.713f, 580.514f, DefaultWidth, DefaultHeight, IdealFisheyeRadius);
    }
};

template<> struct CameraParameters<cammod::Fisheye<double>>
{
    typedef cammod::Fisheye<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 1600;
    static constexpr std::size_t DefaultHeight = 1200;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(315.866, 315.678, 764.713, 580.514, DefaultWidth, DefaultHeight, IdealFisheyeRadius);
    }
};

template<> struct CameraParameters<cammod::Spherical<float>>
{
    typedef cammod::Spherical<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0f, 256.0f, 0.895866f, 1.80258f);
    }
};

template<> struct CameraParameters<cammod::Spherical<double>>
{
    typedef cammod::Spherical<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0, 256.0, 0.895866, 1.80258);
    }
};

template<> struct CameraParameters<cammod::SphericalPovRay<float>>
{
    typedef cammod::SphericalPovRay<float> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0f, 256.0f, 0.895866f, 1.80258f);
    }
};

template<> struct CameraParameters<cammod::SphericalPovRay<double>>
{
    typedef cammod::SphericalPovRay<double> ModelT;
    
    static constexpr std::size_t DefaultWidth = 2048;
    static constexpr std::size_t DefaultHeight = 256;
    
    static void configure(ModelT& camera)
    {
        camera = ModelT(2048.0, 256.0, 0.895866, 1.80258);
    }
};

#endif // CAMERA_PARAMETERS_HPP
