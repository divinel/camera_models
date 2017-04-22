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
 * Tests camera pyramids.
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

#include <CameraModels/CameraModels.hpp>

#include <CameraParameters.hpp>

#include <CameraModels/CameraPyramid.hpp>

template<typename ModelT>
struct IsThisFisheyeModel
{
    static constexpr bool Answer = false;
};

template<> struct IsThisFisheyeModel<cammod::Fisheye<float>> { static constexpr bool Answer = true; };
template<> struct IsThisFisheyeModel<cammod::Fisheye<double>> { static constexpr bool Answer = true; };
template<> struct IsThisFisheyeModel<cammod::FisheyeDistorted<float>> { static constexpr bool Answer = true; };
template<> struct IsThisFisheyeModel<cammod::FisheyeDistorted<double>> { static constexpr bool Answer = true; };

template <typename ModelT>
class PyramidCameraModelTests : public ::testing::Test 
{
public:
    
};

typedef ::testing::Types<
// float
cammod::PinholeDistance<float>,
cammod::PinholeDistanceDistorted<float>,
cammod::PinholeDisparity<float>,
cammod::PinholeDisparityDistorted<float>,
cammod::Generic<float>,
cammod::GenericDistorted<float>,
cammod::Fisheye<float>,
cammod::FisheyeDistorted<float>,
// double
cammod::PinholeDistance<double>,
cammod::PinholeDistanceDistorted<double>,
cammod::PinholeDisparity<double>,
cammod::PinholeDisparityDistorted<double>,
cammod::Generic<double>,
cammod::GenericDistorted<double>,
cammod::Fisheye<double>,
cammod::FisheyeDistorted<double>
> PyramidCameraModelTypes;
TYPED_TEST_CASE(PyramidCameraModelTests, PyramidCameraModelTypes);

TYPED_TEST(PyramidCameraModelTests, TestSimple) 
{
    typedef TypeParam ModelT;
    typedef typename ModelT::Scalar Scalar;
    static constexpr std::size_t PyramidLevels = 2;
    typedef cammod::CameraPyramid<ModelT,PyramidLevels> PyramidT;
    
    ModelT initial_camera;
    CameraParameters<ModelT>::configure(initial_camera);
        
    PyramidT pyr(initial_camera);
    
    typename ModelT::TransformT pose;
    pose.translation() << 10.0 , 20.0, 30.0;
    
    for(std::size_t lvl = 0 ; lvl < PyramidLevels ; ++lvl)
    {
    
        int cnt_good = 0, cnt_bad = 0;
        
        for(unsigned int y = 0 ; y < (unsigned int)pyr[lvl].height() ; ++y)
        {
            for(unsigned int x = 0 ; x < (unsigned int)pyr[lvl].width() ; ++x)
            {
                typename ModelT::PixelT pix((typename ModelT::Scalar)x, (typename ModelT::Scalar)y), pix_out;
                
                // fisheye is not valid outside of the active image area - check pixelValidCircular
                if(!(IsThisFisheyeModel<ModelT>::Answer) || ( pyr[lvl].pixelValidCircular(pix) ) )
                {
                    typename ModelT::PointT pt;
                    const typename ModelT::Scalar distance = 1.5;
                    
                    // lift
                    pt = pyr[lvl].inverseAtDistance(pose, pix, distance); 
                    
                    // project
                    pix_out = pyr[lvl].forward(pose, pt);
                    
                    double err = (pix_out - pix).norm();
                    
                    if(pyr[lvl].pixelValid((Scalar)x,(Scalar)y))
                    {
                        if(err > (Scalar)0.5) // if greater 
                        {
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
        
        EXPECT_EQ(cnt_bad, 0);
    }
}
