/**
 * ****************************************************************************
 * Copyright (c) 2017, Robert Lukierski.
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
 * Radial+Tangential Distortion
 * ****************************************************************************
 */

#ifndef RADIAL_TANGENTIAL_DISTORTION_HPP
#define RADIAL_TANGENTIAL_DISTORTION_HPP

#include <CameraModels/CameraModelUtils.hpp>

namespace cammod
{
/**
 * Full Generic Camera Model, Model Specific Functions.
 * Based on Christopher Mei's model. 
 */
template<typename Derived>
class RadialTangentialDistortion
{
public:
    template<typename T>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE static typename ComplexTypes<T>::PixelT 
        getDistortionVector(const Derived& ccd, const typename ComplexTypes<T>::PixelT& mu)
    {        
        typename ComplexTypes<T>::PixelT ret;
        
        const T x2 = mu(0) * mu(0); // x^2
        const T y2 = mu(1) * mu(1); // y^2
        const T xy = mu(0) * mu(1); // x * y
        const T rho2 = x2 + y2; // rho^2
        const T rho4 = rho2 * rho2; // rho^4
        
        ret(0) = mu(0) * ( ccd.k1() * rho2 + ccd.k2() * rho4 ) + T(2.0) * ccd.p1() * xy + ccd.p2() * ( rho2 + T(2.0) * x2 ); 
        ret(1) = mu(1) * ( ccd.k1() * rho2 + ccd.k2() * rho4 ) + T(2.0) * ccd.p2() * xy + ccd.p1() * ( rho2 + T(2.0) * y2 );
        
        return ret;
    }
    
    template<typename T>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE static typename ComplexTypes<T>::DistortionJacobianT 
        getDistortionPointJacobian(const Derived& ccd, const typename ComplexTypes<T>::PixelT& mu)
    {
        typename ComplexTypes<T>::DistortionJacobianT ret = ComplexTypes<T>::DistortionJacobianT::Zero();
        
        const T x2 = mu(0) * mu(0); // x^2
        const T y2 = mu(1) * mu(1); // y^2
        const T xy = mu(0) * mu(1); // x * y
        const T rho2 = x2 + y2; // rho^2
        const T rho4 = rho2 * rho2; // rho^4
        
        // fwd-x / d{px,py}
        ret(0,0) = T(1.0) + ccd.k1() * rho2 + ccd.k2() * rho4 + T(2.0) * ccd.k1() * x2 + 
                   T(4.0) * ccd.k2() * rho2 * x2 + 
                   T(2.0) * ccd.p1() * mu(1) + 
                   T(6.0) * ccd.p2() * mu(0);
        ret(1,0) = ccd.k1() * T(2.0) * mu(0) * mu(1) + ccd.k2() * T(4.0) * rho2 * mu(0) * mu(1) + 
                   ccd.p1() * T(2.0) * mu(0) + T(2.0) * ccd.p2() * mu(1);
        // fwd-y / d{px,py}
        ret(0,1) = ret(1, 0);
        ret(1,1) = T(1.0) + ccd.k1() * rho2 + ccd.k2() * rho4 + 
                   ccd.k1() * T(2.0) * y2 + ccd.k2() * rho2 * T(4.0) * y2 + 
                   T(6.0) * ccd.p1() * mu(1) + T(2.0) * ccd.p2() * mu(0);
        
        return ret;
    }
    
    template<typename T>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE static typename ComplexTypes<T>::template ForwardParametersJacobianT<4> 
        getDistortionParametersJacobian(const Derived& ccd, const typename ComplexTypes<T>::PixelT& mu)
    {
        typename ComplexTypes<T>::template ForwardParametersJacobianT<4> ret = 
            ComplexTypes<T>::template ForwardParametersJacobianT<4>::Zero();
        
        const T x2 = mu(0) * mu(0); // x^2
        const T y2 = mu(1) * mu(1); // y^2
        const T xy = mu(0) * mu(1); // x * y
        const T rho2 = x2 + y2; // rho^2
        const T rho4 = rho2 * rho2; // rho^4
        
        // fwd-x / d{k1,k2,p1,p2}
        ret(0,0) = mu(0) * rho2;
        ret(0,1) = mu(0) * rho4;
        ret(0,2) = T(2.0) * xy;
        ret(0,3) = y2 + T(3.0) * x2;
        
        // fwd-y / d{k1,k2,p1,p2}
        ret(1,0) = mu(1) * rho2;
        ret(1,1) = mu(1) * rho4;
        ret(1,2) = T(3.0) * y2 + x2;
        ret(1,3) = T(2.0) * xy;
        
        return ret;
    }
};


}

    
#endif // RADIAL_TANGENTIAL_DISTORTION_HP
