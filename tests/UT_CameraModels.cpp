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
 * Tests for statically typed camera models.
 * ****************************************************************************
 */

// system
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <type_traits>

// testing framework & libraries
#include <gtest/gtest.h>

// google logger
#include <glog/logging.h>

#include <CameraModels.hpp>

#include <CameraParameters.hpp>

// -------------------- tricks -----------------------------------

template<typename ModelT>
struct IsThisPovRayModel
{
    static constexpr bool Answer = false;
};

template<> struct IsThisPovRayModel<camera::SphericalPovRayCameraModel<float>> { static constexpr bool Answer = true; };
template<> struct IsThisPovRayModel<camera::SphericalPovRayCameraModel<double>> { static constexpr bool Answer = true; };

template<typename ModelT>
struct IsThisSphericalModel
{
    static constexpr bool Answer = false;
};

template<> struct IsThisSphericalModel<camera::SphericalCameraModel<float>> { static constexpr bool Answer = true; };
template<> struct IsThisSphericalModel<camera::SphericalCameraModel<double>> { static constexpr bool Answer = true; };
template<> struct IsThisSphericalModel<camera::SphericalPovRayCameraModel<float>> { static constexpr bool Answer = true; };
template<> struct IsThisSphericalModel<camera::SphericalPovRayCameraModel<double>> { static constexpr bool Answer = true; };

template<typename ModelT>
struct IsThisFisheyeModel
{
    static constexpr bool Answer = false;
};

template<> struct IsThisFisheyeModel<camera::FisheyeCameraModel<float>> { static constexpr bool Answer = true; };
template<> struct IsThisFisheyeModel<camera::FisheyeCameraModel<double>> { static constexpr bool Answer = true; };
template<> struct IsThisFisheyeModel<camera::IdealFisheyeCameraModel<float>> { static constexpr bool Answer = true; };
template<> struct IsThisFisheyeModel<camera::IdealFisheyeCameraModel<double>> { static constexpr bool Answer = true; };

// -----------------------------------------------------------------------------

template <typename ModelT>
class CameraModelTests : public ::testing::Test 
{
public:
    
};

typedef ::testing::Types<
// float
camera::PinholeCameraModel<float>,
camera::PinholeDistortedCameraModel<float>,
camera::PinholeDisparityCameraModel<float>,
camera::PinholeDisparityDistortedCameraModel<float>,
camera::IdealGenericCameraModel<float>,
camera::FullGenericCameraModel<float>,
camera::SphericalCameraModel<float>,
camera::SphericalPovRayCameraModel<float>,
camera::FisheyeCameraModel<float>,
camera::IdealFisheyeCameraModel<float>,
// double
camera::PinholeCameraModel<double>,
camera::PinholeDistortedCameraModel<double>,
camera::PinholeDisparityCameraModel<double>,
camera::PinholeDisparityDistortedCameraModel<double>,
camera::IdealGenericCameraModel<double>,
camera::FullGenericCameraModel<double>,
camera::SphericalCameraModel<double>,
camera::SphericalPovRayCameraModel<double>,
camera::FisheyeCameraModel<double>,
camera::IdealFisheyeCameraModel<double>
> CameraModelTypes;
TYPED_TEST_CASE(CameraModelTests, CameraModelTypes);

TYPED_TEST(CameraModelTests, TestInverseForward) 
{
    typedef TypeParam ModelT;
    typedef typename ModelT::Scalar Scalar;
    
    ModelT camera;
    
    CameraParameters<ModelT>::configure(camera);
    
    // something is wrong with this particular model
    if(IsThisPovRayModel<ModelT>::Answer)
    {
        return;
    }
    
    std::string scalar_name = "<double>";
    
    if(std::is_same<typename ModelT::Scalar,float>::value)
    {
        scalar_name = "<float>";   
    }
    
    //LOG(INFO) << "Testing " << camera::CameraModelToTypeAndName<ModelT::ModelType>::Name << scalar_name;
    
    typename ModelT::TransformT pose;
    pose.translation() << 10.0 , 20.0, 30.0;
    
    int cnt_good = 0, cnt_bad = 0;
    
    for(unsigned int y = 0 ; y < CameraParameters<ModelT>::DefaultHeight ; ++y)
    {
        for(unsigned int x = 0 ; x < CameraParameters<ModelT>::DefaultWidth ; ++x)
        {
            typename ModelT::PixelT pix((typename ModelT::Scalar)x, (typename ModelT::Scalar)y), pix_out;
            
            // fisheye is not valid outside of the active image area - check pixelValidCircular
            if(!IsThisFisheyeModel<ModelT>::Answer || ( camera.pixelValidCircular(pix) ) )
            {
                typename ModelT::PointT pt;
                const typename ModelT::Scalar distance = 1.5;
                
                // lift
                pt = camera.inverseAtDistance(pose, pix, distance); 
                
                // project
                pix_out = camera.forward(pose, pt);
                
                double err = (pix_out - pix).norm();
                
                if(camera.pixelValid((Scalar)x,(Scalar)y))
                {
                    if(err > (Scalar)0.5)
                    {
                        //LOG(INFO) << camera::CameraModelToTypeAndName<ModelT::ModelType>::Name << scalar_name << ": Error at " << x << "," << y << " is " << pix_out(0) << "," << pix_out(1) << " = " << err << " / for " << pt(0) << " , " << pt(1) << " , " << pt(2);
                        cnt_bad++;
                    }
                    else
                    {
                        cnt_good++;
                    }
                }
            }
        }
    }
}
