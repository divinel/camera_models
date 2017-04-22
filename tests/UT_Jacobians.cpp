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
 * Tests for hardcoded Jacobians.
 * ****************************************************************************
 */

// system
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <type_traits>
#include <chrono>

// testing framework & libraries
#include <gtest/gtest.h>

// google logger
#include <glog/logging.h>

#include <CameraModels/CameraModels.hpp>

#include <CameraParameters.hpp>

#include <unsupported/Eigen/AutoDiff>

template <typename ModelT>
class JacobianCameraModelTests : public ::testing::Test 
{
public:
    
};

typedef ::testing::Types<
// double
cammod::PinholeDistance<double>,
cammod::PinholeDistanceDistorted<double>,
cammod::PinholeDisparity<double>,
cammod::PinholeDisparityDistorted<double>,
cammod::Generic<double>,
cammod::GenericDistorted<double>,
cammod::Fisheye<double>
> JacobianCameraModelTypes;

TYPED_TEST_CASE(JacobianCameraModelTests, JacobianCameraModelTypes);

TYPED_TEST(JacobianCameraModelTests, TestForwardPoint) 
{
    typedef TypeParam ModelT;
    typedef typename ModelT::Scalar Scalar;
    typedef typename cammod::ComplexTypes<Scalar>::PointT PointT;
    typedef typename cammod::ComplexTypes<Scalar>::ForwardPointJacobianT JacobianT;
    
    ModelT tmp_camera;
    CameraParameters<ModelT>::configure(tmp_camera);
    
    typedef Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,3,1>> JetT;
    typedef typename cammod::ComplexTypes<JetT>::PixelT JetPixelT;
    typedef typename cammod::ComplexTypes<JetT>::PointT JetPointT;
    
    const PointT ptin(1.0,2.0,3.0);
    const JetPointT jptin(JetT(ptin(0),3,0),JetT(ptin(1),3,1),JetT(ptin(2),3,2));

    const JetPixelT jpix = tmp_camera.template forward<JetT>(jptin);
    const JacobianT jhard = tmp_camera.template forwardPointJacobian(ptin);
    JacobianT jauto;
    jauto.template block<1,3>(0,0) = jpix(0).derivatives();
    jauto.template block<1,3>(1,0) = jpix(1).derivatives();
    
    //std::cerr << "Point Jacobian Auto: " << std::endl << jauto << std::endl;
    //std::cerr << "Point Jacobian Hardcoded: " << std::endl << jhard << std::endl;
    //std::cerr << "Jacobian Error: " << (jauto - jhard).norm() << std::endl;
    
    for(std::size_t c = 0 ; c < 3 ; ++c)
    {
        for(std::size_t r = 0 ; r < 2 ; ++r)
        {
            EXPECT_FLOAT_EQ(jauto(r,c), jhard(r,c)) << "Point Jacobian not equal at (" << r << "," << c << ")";
        }
    }
}

TYPED_TEST(JacobianCameraModelTests, TestForwardParameters) 
{
    typedef TypeParam ModelT;
    typedef typename ModelT::Scalar Scalar;
    typedef typename cammod::ComplexTypes<Scalar>::PointT PointT;
    typedef typename ModelT::template ForwardParametersJacobianT<Scalar> JacobianT;
    typedef Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,ModelT::ParametersToOptimize,1>> JetT;
    typedef typename cammod::ComplexTypes<JetT>::PixelT JetPixelT;
    typedef typename cammod::ComplexTypes<JetT>::PointT JetPointT;
    typedef typename ModelT::template GetOtherType<JetT> JetModelT;
    
    ModelT tmp_camera;
    CameraParameters<ModelT>::configure(tmp_camera);
    JetModelT jtmp_camera;
    
    for(std::size_t i = 0 ; i < ModelT::ParametersToOptimize ; ++i)
    {
        jtmp_camera.data()[i] = JetT(tmp_camera.data()[i],ModelT::ParametersToOptimize,i);
    }
    
    const PointT ptin(1.0,2.0,3.0);
    const JetPointT jptin(JetT(ptin(0)),JetT(ptin(1)),JetT(ptin(2)));

    const JetPixelT jpix = jtmp_camera.template forward<JetT>(jptin);
    const JacobianT jhard = tmp_camera.template forwardParametersJacobian(ptin);
    JacobianT jauto;
    
    jauto.template block<1,ModelT::ParametersToOptimize>(0,0) = jpix(0).derivatives();
    jauto.template block<1,ModelT::ParametersToOptimize>(1,0) = jpix(1).derivatives();
    
    //std::cerr << "Parameters Jacobian Auto: " << std::endl << jauto << std::endl;
    //std::cerr << "Parameters Jacobian Hardcoded: " << std::endl << jhard << std::endl;
    //std::cerr << "Jacobian Error: " << (jauto - jhard).norm() << std::endl;
    
    for(std::size_t c = 0 ; c < ModelT::ParametersToOptimize ; ++c)
    {
        for(std::size_t r = 0 ; r < 2 ; ++r)
        {
            EXPECT_FLOAT_EQ(jauto(r,c), jhard(r,c)) << "Parameters  Jacobian not equal at (" << r << "," << c << ")";
        }
    }
}
