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
 * Distorted Gayer-Baretto-Mei generic camera model.
 * ****************************************************************************
 */

#ifndef GENERIC_DISTORTED_CAMERA_MODEL_HPP
#define GENERIC_DISTORTED_CAMERA_MODEL_HPP

#include <CameraModels/CameraModelUtils.hpp>
#include <CameraModels/Generic.hpp>

// fwd
namespace cammod 
{
template<typename _Scalar, int _Options = 0> class GenericDistorted;

namespace internal
{
static constexpr unsigned int GenericDistortedParameterCount = 14;
}
}

// Eigen Traits 
namespace Eigen 
{
    namespace internal 
    {
        template<typename _Scalar, int _Options>
        struct traits<cammod::GenericDistorted<_Scalar,_Options> > 
        {
            typedef _Scalar Scalar;
            typedef Matrix<Scalar,cammod::internal::GenericDistortedParameterCount,1> ComplexType;
            static constexpr bool HasForwardPointJacobian = true;
            static constexpr bool HasForwardParametersJacobian = true;
            static constexpr bool HasInversePointJacobian = false;
            static constexpr bool HasInverseParametersJacobian = false;
            static constexpr unsigned int NumParameters = cammod::internal::GenericDistortedParameterCount;
            static constexpr unsigned int ParametersToOptimize = NumParameters - 4;
            static constexpr bool CalibrationSupported = true;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<cammod::GenericDistorted<_Scalar>, _Options> > : traits<cammod::GenericDistorted<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<Matrix<Scalar,cammod::internal::GenericDistortedParameterCount,1>,_Options> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<const cammod::GenericDistorted<_Scalar>, _Options> > : traits<const cammod::GenericDistorted<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<const Matrix<Scalar,cammod::internal::GenericDistortedParameterCount,1>,_Options> ComplexType;
        };    
    }
}

namespace cammod
{
/**
 * Full Generic Camera Model, Model Specific Functions.
 * Based on Christopher Mei's model. 
 */
template<typename Derived>
class GenericDistortedBase : public CameraFunctions<Derived>, public ComplexTypes<typename Eigen::internal::traits<Derived>::Scalar>
{
    typedef CameraFunctions<Derived> FunctionsBase;
public:
    // various helpers and types
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
    typedef typename Eigen::internal::traits<Derived>::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Derived>::ComplexType& ConstComplexReference;
    static constexpr CameraModelType ModelType = CameraModelType::GenericDistorted;
    
    using FunctionsBase::forward;
    using FunctionsBase::inverse;
    using FunctionsBase::inverseAtDistance;
    using FunctionsBase::forwardPointJacobian;
    using FunctionsBase::forwardParametersJacobian;
    using FunctionsBase::twoFrameProject;
    using FunctionsBase::worldToCamera;
    using FunctionsBase::cameraToWorld;
    using FunctionsBase::pointValid;
    using FunctionsBase::pixelValid;
    using FunctionsBase::pixelValidSquare;
    using FunctionsBase::pixelValidCircular;
    using FunctionsBase::resizeViewport;
    
    template<typename NewScalarType>
    EIGEN_DEVICE_FUNC inline GenericDistorted<NewScalarType> cast() const { return GenericDistorted<NewScalarType>(access().template cast<NewScalarType>()); }
    
    template<typename OtherDerived> 
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE GenericDistortedBase<Derived>& operator=(const GenericDistortedBase<OtherDerived> & other) { access_nonconst() = other.access(); return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const  { return static_cast<const Derived*>(this)->access(); }
private:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst()  { return static_cast<Derived*>(this)->access_nonconst(); }
public:
    static constexpr unsigned int NumParameters = Eigen::internal::traits<Derived>::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Eigen::internal::traits<Derived>::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Eigen::internal::traits<Derived>::CalibrationSupported;
    static constexpr bool HasForwardPointJacobian = Eigen::internal::traits<Derived>::HasForwardPointJacobian;
    static constexpr bool HasForwardParametersJacobian = Eigen::internal::traits<Derived>::HasForwardParametersJacobian;
    static constexpr bool HasInversePointJacobian = Eigen::internal::traits<Derived>::HasInversePointJacobian;
    static constexpr bool HasInverseParametersJacobian = Eigen::internal::traits<Derived>::HasInverseParametersJacobian;
    
    template<typename T = Scalar>
    using ForwardParametersJacobianT = typename ComplexTypes<T>::template ForwardParametersJacobianT<ParametersToOptimize>;
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverse(const Derived& ccd, T x, T y) 
    {
        typename ComplexTypes<T>::PointT ret;
        
        // inverse intrinsics pixel -> image plane
        typename ComplexTypes<T>::PixelT pix_udist, pix_dist;
        pix_udist(0) = (ccd.fx() * ccd.skew() * ccd.v0() - ccd.fy() * ccd.u0())/(ccd.fx() * ccd.fy()) - (y * ccd.skew())/ccd.fy() + x / ccd.fx();
        pix_udist(1) = (y - ccd.v0()) / ccd.fy();
        
        // inverse distortion - distorted image plane
        pix_dist = pix_udist - getDistortionVector<T>(ccd, pix_udist - getDistortionVector<T>(ccd, pix_udist - getDistortionVector<T>(ccd, pix_udist - getDistortionVector<T>(ccd, pix_udist - getDistortionVector<T>(ccd, pix_udist))))); 
        
        // inverse perspective - pixel to point
        const T term = (ccd.epsilon() + sqrt(T(1.0f) + (T(1.0f) - ccd.epsilon() * ccd.epsilon()) * (pix_dist(0) * pix_dist(0) + pix_dist(1) * pix_dist(1)))) 
                       / (pix_dist(0) * pix_dist(0) + pix_dist(1) * pix_dist(1) + T(1.0f));
        
        ret(0) = term * pix_dist(0);
        ret(1) = term * pix_dist(1); 
        ret(2) = term - ccd.epsilon();
        
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt) 
    {
        typename ComplexTypes<T>::PixelT ret, p;
        
        // unit vector
        const typename ComplexTypes<T>::PointT unit_pt = tmp_pt.normalized(); 
        
        // perspective
        p = unit_pt.template topRows<2>() / (unit_pt(2) + ccd.epsilon());

        // distortions
        p += getDistortionVector<T>(ccd, p);
        
        // intrinsics
        ret(0) = ccd.fx() * p(0) + (ccd.fx() * ccd.skew()) * p(1) + ccd.u0();
        ret(1) = ccd.fy() * p(1) + ccd.v0();
        
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::ForwardPointJacobianT forwardPointJacobian(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt)
    {
        using Eigen::numext::sqrt;
        
        typename ComplexTypes<T>::ForwardPointJacobianT perspective_jacobian = ComplexTypes<T>::ForwardPointJacobianT::Zero();
        const T x2 = tmp_pt(0) * tmp_pt(0);
        const T y2 = tmp_pt(1) * tmp_pt(1);
        const T z2 = tmp_pt(2) * tmp_pt(2);
        const T xy = tmp_pt(0) * tmp_pt(1);
        const T r = tmp_pt.norm();
        const T r2 = tmp_pt.squaredNorm();
        const T denom = (tmp_pt(2)/r + ccd.epsilon());
        const T denom1 = sqrt(r2*r2*r2) * denom;
        const T denom2 = r2 * r2 * denom * denom;
        const T denom3 = r * denom * denom;
        const T term1 = z2 / sqrt(r2*r2*r2);
        
        // perspective-x / d{x,y,z}
        perspective_jacobian(0,0) = T(1.0) / (r * denom) - x2 /  denom1 + (x2 * tmp_pt(2)) / denom2;
        perspective_jacobian(0,1) = (xy * tmp_pt(2)) / denom2 - xy / denom1;
        perspective_jacobian(0,2) =-(tmp_pt(0) * tmp_pt(2)) / denom1 - ( tmp_pt(0) * (T(1.0) / r - term1 ) ) / denom3;
        
        // perspective-y / d{x,y,z}
        perspective_jacobian(1,0) = (xy * tmp_pt(2)) / denom2 - xy / denom1;
        perspective_jacobian(1,1) = T(1.0) / (r * denom) - y2 /  denom1 + (y2 * tmp_pt(2)) / denom2;
        perspective_jacobian(1,2) =-(tmp_pt(1) * tmp_pt(2)) / denom1 - ( tmp_pt(1) * (T(1.0) / r - term1 ) ) / denom3;
        
        // distortion / d{px,py}
        const typename ComplexTypes<T>::PointT unit_pt = tmp_pt.normalized();        
        const typename ComplexTypes<T>::PixelT undist_pt(unit_pt(0)/(unit_pt(2) + ccd.epsilon()),unit_pt(1)/(unit_pt(2) + ccd.epsilon()));
        const typename ComplexTypes<T>::DistortionJacobianT distortion_jacobian = getDistortionPointJacobian<T>(ccd, undist_pt);
        
        // intrinsics / d{px,py}
        typename ComplexTypes<T>::template ForwardParametersJacobianT<2> intrinsics_jacobian = ComplexTypes<T>::template ForwardParametersJacobianT<2>::Zero();
        intrinsics_jacobian(0,0) = ccd.fx();
        intrinsics_jacobian(0,1) = ccd.fx() * ccd.skew();
        intrinsics_jacobian(1,1) = ccd.fy();
        
        return intrinsics_jacobian * distortion_jacobian * perspective_jacobian;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ForwardParametersJacobianT<T> forwardParametersJacobian(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt)
    {
        ForwardParametersJacobianT<T> ret = ForwardParametersJacobianT<T>::Zero();
        
        const typename ComplexTypes<T>::PointT unit_pt = tmp_pt.normalized();        
        const typename ComplexTypes<T>::PixelT undist_pt(unit_pt(0)/(unit_pt(2) + ccd.epsilon()),unit_pt(1)/(unit_pt(2) + ccd.epsilon()));
        const typename ComplexTypes<T>::PixelT dist_pt = undist_pt + getDistortionVector<T>(ccd, undist_pt);
     
        // dforward / d{epsilon}
        typename ComplexTypes<T>::template ForwardParametersJacobianT<1> epsilon_jacobian = ComplexTypes<T>::template ForwardParametersJacobianT<1>::Zero();
        const T r = tmp_pt.norm();
        const T denom = (tmp_pt(2) + ccd.epsilon() * r) * (tmp_pt(2) + ccd.epsilon() * r);
        epsilon_jacobian(0,0) = - (r * tmp_pt(0)) / denom;
        epsilon_jacobian(1,0) = - (r * tmp_pt(1)) / denom;
        
        // dforward-distortions / d{k1,k2,p1,p2}
        const typename ComplexTypes<T>::template ForwardParametersJacobianT<4> distortions_jacobian = getDistortionParametersJacobian<T>(ccd,undist_pt);
        
        // dforward-distortions / d{px,py}
        const typename ComplexTypes<T>::DistortionJacobianT distortions_to_chain = getDistortionPointJacobian<T>(ccd,undist_pt);
        
        // dforward-intrinsics / d{fx,fy,u0,v0,skew}
        typename ComplexTypes<T>::template ForwardParametersJacobianT<5> intrinsics_jacobian = ComplexTypes<T>::template ForwardParametersJacobianT<5>::Zero();
        intrinsics_jacobian(0,0) = dist_pt(1) * ccd.skew() + dist_pt(0); 
        intrinsics_jacobian(0,2) = T(1.0);
        intrinsics_jacobian(0,4) = dist_pt(1) * ccd.fx();  // skew
        intrinsics_jacobian(1,1) = dist_pt(1);
        intrinsics_jacobian(1,3) = T(1.0);
        
        // dforward-intrinsics / d{px,py}
        typename ComplexTypes<T>::template ForwardParametersJacobianT<2> intrinsics_to_chain = ComplexTypes<T>::template ForwardParametersJacobianT<2>::Zero();
        intrinsics_to_chain(0,0) = ccd.fx();
        intrinsics_to_chain(0,1) = ccd.fx() * ccd.skew();
        intrinsics_to_chain(1,1) = ccd.fy();

        ret.template block<2,4>(0,0) = intrinsics_jacobian.template block<2,4>(0,0); 
        ret.template block<2,1>(0,4) = (intrinsics_to_chain * (distortions_to_chain * epsilon_jacobian)); 
        ret.template block<2,4>(0,5) = intrinsics_to_chain * distortions_jacobian; 
        ret.template block<2,1>(0,9) = intrinsics_jacobian.template block<2,1>(0,4);
        
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pointValid(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt)
    {
        using Eigen::numext::abs;
        // don't divide by zero
        const typename ComplexTypes<T>::PointT unit_pt = tmp_pt.normalized();
        return abs(unit_pt(2)) > Eigen::NumTraits<T>::epsilon();
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidSquare(const Derived& ccd, T x, T y) 
    {
        if((x >= T(0.0)) && (x < ccd.width()) && (y >= T(0.0)) && (y < ccd.height()))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidCircular(const Derived& ccd, T x, T y) 
    {
        if(ccd.r1() <= 0.0) // no r1
        {
            const T x2 = (x - ccd.u0()) * (x - ccd.u0());
            const T y2 = (y - ccd.v0()) * (y - ccd.v0());
            if(x2 + y2 < ccd.r2() * ccd.r2())
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            const T x2 = (x - ccd.u0()) * (x - ccd.u0());
            const T y2 = (y - ccd.v0()) * (y - ccd.v0());
            const T r12 = ccd.r1() * ccd.r1();
            const T r22 = ccd.r2() * ccd.r2();
            
            if(((x2 + y2) > r12) && ((x2 + y2) < r22))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void resizeViewport(Derived& ccd, const T& new_width, const T& new_height)
    {
        using Eigen::numext::mini;
        
        const T x_ratio = new_width / ccd.width();
        const T y_ratio = new_height / ccd.height();
        const T r_ratio = mini(x_ratio, y_ratio);
        
        ccd.fx(ccd.fx() * x_ratio);
        ccd.fy(ccd.fy() * y_ratio);
        ccd.u0(ccd.u0() * x_ratio);
        ccd.v0(ccd.v0() * y_ratio);
        ccd.r1(ccd.r1() * r_ratio);
        ccd.r2(ccd.r2() * r_ratio);
        ccd.width(new_width);
        ccd.height(new_height);
    }
    
    // access to camera model parameters
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& fx() const { return data()[0]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& fy() const { return data()[1]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& u0() const { return data()[2]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& v0() const { return data()[3]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& epsilon() const { return data()[4]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k1() const { return data()[5]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k2() const { return data()[6]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& p1() const { return data()[7]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& p2() const { return data()[8]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& skew() const { return data()[9]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& width() const { return data()[10]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& height() const { return data()[11]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& r1() const { return data()[12]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& r2() const { return data()[13]; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fx(const Scalar& v) { data()[0] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fy(const Scalar& v) { data()[1] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void u0(const Scalar& v) { data()[2] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void v0(const Scalar& v) { data()[3] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void epsilon(const Scalar& v) { data()[4] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k1(const Scalar& v) { data()[5] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k2(const Scalar& v) { data()[6] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void p1(const Scalar& v) { data()[7] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void p2(const Scalar& v) { data()[8] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void skew(const Scalar& v) { data()[9] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void width(const Scalar& v) { data()[10] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void height(const Scalar& v) { data()[11] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void r1(const Scalar& v) { data()[12] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void r2(const Scalar& v) { data()[13] = v; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar* data() { return access_nonconst().data(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar* data() const { return access().data(); }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Generic<Scalar> getIdeal() const
    {
        return Generic<Scalar>(fx(),fy(),u0(),v0(),epsilon(),width(),height(),r1(),r2());
    }
    
#ifdef CAMERA_MODELS_HAVE_SERIALIZER
    template <typename Archive>
    void serialize(Archive & archive, std::uint32_t const version)
    { 
        CAMERA_MODELS_SERIALIZE(archive,"fx",data()[0]);
        CAMERA_MODELS_SERIALIZE(archive,"fy",data()[1]);
        CAMERA_MODELS_SERIALIZE(archive,"u0",data()[2]);
        CAMERA_MODELS_SERIALIZE(archive,"v0",data()[3]);
        CAMERA_MODELS_SERIALIZE(archive,"epsilon",data()[4]);
        CAMERA_MODELS_SERIALIZE(archive,"k1",data()[5]);
        CAMERA_MODELS_SERIALIZE(archive,"k2",data()[6]);
        CAMERA_MODELS_SERIALIZE(archive,"p1",data()[7]);
        CAMERA_MODELS_SERIALIZE(archive,"p2",data()[8]);
        CAMERA_MODELS_SERIALIZE(archive,"skew",data()[9]);
        CAMERA_MODELS_SERIALIZE(archive,"width",data()[10]);
        CAMERA_MODELS_SERIALIZE(archive,"height",data()[11]);
        CAMERA_MODELS_SERIALIZE(archive,"r1",data()[12]);
        CAMERA_MODELS_SERIALIZE(archive,"r2",data()[13]);
    }
#endif // CAMERA_MODELS_HAVE_SERIALIZER
    
private:
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC static inline typename ComplexTypes<T>::PixelT getDistortionVector(const Derived& ccd, const typename ComplexTypes<T>::PixelT& mu)
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
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC static inline typename ComplexTypes<T>::DistortionJacobianT getDistortionPointJacobian(const Derived& ccd, const typename ComplexTypes<T>::PixelT& mu)
    {
        typename ComplexTypes<T>::DistortionJacobianT ret = ComplexTypes<T>::DistortionJacobianT::Zero();
        
        const T x2 = mu(0) * mu(0); // x^2
        const T y2 = mu(1) * mu(1); // y^2
        const T xy = mu(0) * mu(1); // x * y
        const T rho2 = x2 + y2; // rho^2
        const T rho4 = rho2 * rho2; // rho^4
        
        // fwd-x / d{px,py}
        ret(0,0) = T(1.0) + ccd.k1() * rho2 + ccd.k2() * rho4 + T(2.0) * ccd.k1() * x2 + T(4.0) * ccd.k2() * rho2 * x2 + T(2.0) * ccd.p1() * mu(1) + T(6.0) * ccd.p2() * mu(0);
        ret(1,0) = ccd.k1() * T(2.0) * mu(0) * mu(1) + ccd.k2() * T(4.0) * rho2 * mu(0) * mu(1) + ccd.p1() * T(2.0) * mu(0) + T(2.0) * ccd.p2() * mu(1);
        // fwd-y / d{px,py}
        ret(0,1) = ret(1, 0);
        ret(1,1) = T(1.0) + ccd.k1() * rho2 + ccd.k2() * rho4 + ccd.k1() * T(2.0) * y2 + ccd.k2() * rho2 * T(4.0) * y2 + T(6.0) * ccd.p1() * mu(1) + T(2.0) * ccd.p2() * mu(0);
        
        return ret;
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC static inline typename ComplexTypes<T>::template ForwardParametersJacobianT<4> getDistortionParametersJacobian(const Derived& ccd, const typename ComplexTypes<T>::PixelT& mu)
    {
        typename ComplexTypes<T>::template ForwardParametersJacobianT<4> ret = ComplexTypes<T>::template ForwardParametersJacobianT<4>::Zero();
        
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

/**
 * Full Generic Camera Model, Eigen storage.
 */
template<typename _Scalar, int _Options>
class GenericDistorted : public GenericDistortedBase<GenericDistorted<_Scalar, _Options>>
{
    typedef GenericDistortedBase<GenericDistorted<_Scalar, _Options>> Base;
public:
    typedef typename Eigen::internal::traits<GenericDistorted<_Scalar,_Options> >::Scalar Scalar;
    typedef typename Eigen::internal::traits<GenericDistorted<_Scalar,_Options> >::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<GenericDistorted<_Scalar,_Options> >::ComplexType& ConstComplexReference;
    
    template<class OtherT> using GetOtherType = GenericDistorted<OtherT, _Options>;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class GenericDistortedBase<GenericDistorted<_Scalar, _Options>>;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    EIGEN_DEVICE_FUNC inline GenericDistorted(Scalar fx, Scalar fy, Scalar u0, Scalar v0, Scalar epsilon = 0.0, Scalar k1 = 0.0, Scalar k2 = 0.0, Scalar p1 = 0.0, Scalar p2 = 0.0, Scalar skew = 0.0, Scalar w = 0.0, Scalar h = 0.0, Scalar r1 = 0.0, Scalar r2 = 0.0)
    {
        access_nonconst() << fx , fy , u0 , v0 , epsilon , k1 , k2 , p1 , p2 , skew , w , h , r1 , r2;
    }
    
    EIGEN_DEVICE_FUNC inline GenericDistorted() : parameters(Eigen::Matrix<Scalar,NumParameters,1>::Zero()) { }
    EIGEN_DEVICE_FUNC inline GenericDistorted(const typename Eigen::internal::traits<GenericDistorted<_Scalar,_Options> >::ComplexType& vec) : parameters(vec) { }
    
    EIGEN_DEVICE_FUNC inline GenericDistorted& operator=(const typename Eigen::internal::traits<GenericDistorted<_Scalar,_Options> >::ComplexType& vec) { access_nonconst() = vec; return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    Eigen::Matrix<Scalar,NumParameters,1> parameters;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const GenericDistorted<T>& p)
{
    os << "FullGeneric(fx = " << p.fx() << ", fy = " << p.fy() << ", u0 = " << p.u0() << ", v0 = " << p.v0() << ", eps = " << p.epsilon() << ", k1 = " << p.k1() << ", k2 = " << p.k2() << ", p1 = " << p.p1() << ", p2 = " << p.p2() << ", s = " << p.skew() << ", " << p.width() << " x " << p.height() << ", r1 = " << p.r1() << ", r2 = " << p.r2() << ")";
    return os;
}

}

namespace Eigen 
{
/**
 * Full Generic Camera Model, Eigen Map.
 */
template<typename _Scalar, int _Options>
class Map<cammod::GenericDistorted<_Scalar>, _Options> : public cammod::GenericDistortedBase<Map<cammod::GenericDistorted<_Scalar>, _Options>>
{
    typedef cammod::GenericDistortedBase<Map<cammod::GenericDistorted<_Scalar>, _Options>> Base;
    
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef typename internal::traits<Map>::ComplexType& ComplexReference;
    typedef const typename internal::traits<Map>::ComplexType& ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class cammod::GenericDistortedBase<Map<cammod::GenericDistorted<_Scalar>, _Options>>;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(Scalar* coeffs) : parameters(coeffs)  { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    
    Map<Matrix<Scalar,NumParameters,1>,_Options> parameters;
};

/**
 * Full Generic Camera Model, Eigen Map const.
 */
template<typename _Scalar, int _Options>
class Map<const cammod::GenericDistorted<_Scalar>, _Options> : public cammod::GenericDistortedBase<Map<const cammod::GenericDistorted<_Scalar>, _Options>>
{
    typedef cammod::GenericDistortedBase<Map<const cammod::GenericDistorted<_Scalar>, _Options>> Base;
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef const typename internal::traits<Map>::ComplexType & ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(const Scalar* coeffs) : parameters(coeffs) { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const  { return parameters; }
protected:
    const Map<const Matrix<Scalar,NumParameters,1>,_Options> parameters;
};
}
    
#endif // GENERIC_DISTORTED_CAMERA_MODEL_HPP
