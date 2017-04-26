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
 * Distorted Pinhole Camera model with inverse in distance.
 * ****************************************************************************
 */

#ifndef PINHOLE_DISTANCE_DISTORTED_CAMERA_MODEL_HPP
#define PINHOLE_DISTANCE_DISTORTED_CAMERA_MODEL_HPP

#include <CameraModels/CameraModelUtils.hpp>
#include <CameraModels/PinholeDistance.hpp>
#include <CameraModels/RadialTangentialDistortion.hpp>

// fwd
namespace cammod 
{
template<typename _Scalar, int _Options = 0> class PinholeDistanceDistorted;

namespace internal
{
static constexpr unsigned int PinholeDistanceDistortedParameterCount = 11;
}
}

// Eigen Traits, but also some traits for our use
namespace Eigen 
{
    namespace internal 
    {
        template<typename _Scalar, int _Options>
        struct traits<cammod::PinholeDistanceDistorted<_Scalar,_Options> > 
        {
            typedef _Scalar Scalar;
            typedef Matrix<Scalar,cammod::internal::PinholeDistanceDistortedParameterCount,1> ComplexType;
            static constexpr bool HasForwardPointJacobian = true;
            static constexpr bool HasForwardParametersJacobian = true;
            static constexpr bool HasInversePointJacobian = false;
            static constexpr bool HasInverseParametersJacobian = false;
            static constexpr unsigned int NumParameters = cammod::internal::PinholeDistanceDistortedParameterCount;
            static constexpr unsigned int ParametersToOptimize = NumParameters - 2;
            static constexpr bool CalibrationSupported = true;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<cammod::PinholeDistanceDistorted<_Scalar>, _Options> >
            : traits<cammod::PinholeDistanceDistorted<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<Matrix<Scalar,cammod::internal::PinholeDistanceDistortedParameterCount,1>,_Options> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<const cammod::PinholeDistanceDistorted<_Scalar>, _Options> >
            : traits<const cammod::PinholeDistanceDistorted<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<const Matrix<Scalar,cammod::internal::PinholeDistanceDistortedParameterCount,1>,_Options> ComplexType;
        };    
    }
}

namespace cammod
{
/**
 * PinholeDistanceDistorted, Model Specific Functions.
 */
template<typename Derived>
class PinholeDistanceDistortedBase : public CameraFunctions<Derived>, 
                                     public ComplexTypes<typename Eigen::internal::traits<Derived>::Scalar>
{
    typedef CameraFunctions<Derived> FunctionsBase;
    typedef RadialTangentialDistortion<Derived> DistortionT;
public:
    // various helpers and types
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
    typedef typename Eigen::internal::traits<Derived>::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Derived>::ComplexType& ConstComplexReference;
    static constexpr CameraModelType ModelType = CameraModelType::PinholeDistanceDistorted;
    
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
    EIGEN_DEVICE_FUNC inline PinholeDistanceDistorted<NewScalarType> cast() const 
    { 
        return PinholeDistanceDistorted<NewScalarType>(access().template cast<NewScalarType>()); 
    }
    
    template<typename OtherDerived> 
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE PinholeDistanceDistortedBase<Derived>& 
        operator=(const PinholeDistanceDistortedBase<OtherDerived> & other) 
    { 
        access_nonconst() = other.access(); 
        return *this; 
    }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const 
    { 
        return static_cast<const Derived*>(this)->access(); 
    }
private:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() 
    { 
        return static_cast<Derived*>(this)->access_nonconst();       
    }
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
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT 
        inverse(const Derived& ccd, T x, T y) 
    {
        using Eigen::numext::sqrt;
        
        typename ComplexTypes<T>::PointT ret;
        
        // inverse intrinsics pixel -> image plane
        typename ComplexTypes<T>::PixelT pix_udist, pix_dist;
        pix_udist(0) = (ccd.fx() * ccd.skew() * ccd.v0() - ccd.fy() * ccd.u0()) / 
                       (ccd.fx() * ccd.fy()) - (y * ccd.skew())/ccd.fy() + x / ccd.fx();
        pix_udist(1) = (y - ccd.v0()) / ccd.fy();
        
        // inverse distortion - distorted image plane
        pix_dist = pix_udist - DistortionT::template getDistortionVector<T>(ccd, pix_udist - DistortionT::template getDistortionVector<T>(ccd, 
                                                                                 pix_udist - DistortionT::template getDistortionVector<T>(ccd, 
                                                                                 pix_udist - DistortionT::template getDistortionVector<T>(ccd, pix_udist)))); 
      
        const T xx = pix_dist(0);
        const T yy = pix_dist(1);
        const T z =  T(1.0) / sqrt(xx*xx + yy*yy + T(1.0) );
        
        ret(0) = xx*z;
        ret(1) = yy*z;
        ret(2) = z;

        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT 
        forward(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt) 
    {
        typename ComplexTypes<T>::PixelT ret, p;
        
        // perspective        
        p(0) = tmp_pt(0) / tmp_pt(2);
        p(1) = tmp_pt(1) / tmp_pt(2);

        // distortions
        p += DistortionT::template getDistortionVector<T>(ccd, p);
        
        // intrinsics
        ret(0) = ccd.fx() * p(0) + (ccd.fx() * ccd.skew()) * p(1) + ccd.u0();
        ret(1) = ccd.fy() * p(1) + ccd.v0();
        
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::ForwardPointJacobianT forwardPointJacobian(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt)
    {
        typename ComplexTypes<T>::ForwardPointJacobianT perspective_jacobian = ComplexTypes<T>::ForwardPointJacobianT::Zero();
        typename ComplexTypes<T>::template ForwardParametersJacobianT<2> intrinsics_jacobian = 
            ComplexTypes<T>::template ForwardParametersJacobianT<2>::Zero();
        
        const T invz = T(1.0) / tmp_pt(2);
        const T invz2 = invz * invz;
        
        // perspective / d{x,y,z}
        perspective_jacobian(0,0) = invz;
        perspective_jacobian(0,2) = -tmp_pt(0) * invz2;
        perspective_jacobian(1,1) = invz;
        perspective_jacobian(1,2) = -tmp_pt(1) * invz2;
        
        // distortion / d{px,py}
        const typename ComplexTypes<T>::PixelT undist_pt(tmp_pt(0)/tmp_pt(2),tmp_pt(1)/tmp_pt(2));
        const typename ComplexTypes<T>::DistortionJacobianT distortion_jacobian = 
            DistortionT::template getDistortionPointJacobian<T>(ccd, undist_pt);
        
        // intrinsics / d{px,py}
        intrinsics_jacobian(0,0) = ccd.fx();
        intrinsics_jacobian(0,1) = ccd.fx() * ccd.skew();
        intrinsics_jacobian(1,1) = ccd.fy();
        
        return intrinsics_jacobian * distortion_jacobian * perspective_jacobian;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ForwardParametersJacobianT<T> 
        forwardParametersJacobian(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt)
    {
        ForwardParametersJacobianT<T> ret = ForwardParametersJacobianT<T>::Zero();
        
        // intrinsics / d{px,py}
        typename ComplexTypes<T>::template ForwardParametersJacobianT<2> intrinsics_jacobian = 
            ComplexTypes<T>::template ForwardParametersJacobianT<2>::Zero();
        intrinsics_jacobian(0,0) = ccd.fx();
        intrinsics_jacobian(0,1) = ccd.fx() * ccd.skew();
        intrinsics_jacobian(1,1) = ccd.fy();
        
        // distortions / d{k1,k2,p1,p2}
        const typename ComplexTypes<T>::PixelT undist_pt(tmp_pt(0)/tmp_pt(2),tmp_pt(1)/tmp_pt(2));
        const typename ComplexTypes<T>::PixelT dist_pt = undist_pt + DistortionT::template getDistortionVector<T>(ccd, undist_pt);
        const typename ComplexTypes<T>::template ForwardParametersJacobianT<4> dist_jac = 
            DistortionT::template getDistortionParametersJacobian<T>(ccd,undist_pt);
        
        // int-fwd-x / d{k1,k2,p1,p2,skew}
        ret(0,0) = dist_pt(1) * ccd.skew() + dist_pt(0); 
        ret(0,2) = T(1.0);
        ret(0,8) = dist_pt(1) * ccd.fx(); 
        
        // int-fwd-y / d{k1,k2,p1,p2,skew}
        ret(1,1) = dist_pt(1);
        ret(1,3) = T(1.0);
        
        // plug distortions with chain rule
        ret.template block<2,4>(0,4) = (intrinsics_jacobian * dist_jac);
        
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pointValid(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt)
    {
        using Eigen::numext::abs;
        // don't divide by zero
        return abs(tmp_pt(2)) > Eigen::NumTraits<T>::epsilon();
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
        return true;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void resizeViewport(Derived& ccd, const T& new_width, const T& new_height)
    {
        const T x_ratio = new_width / ccd.width();
        const T y_ratio = new_height / ccd.height();
        
        ccd.fx(ccd.fx() * x_ratio);
        ccd.fy(ccd.fy() * y_ratio);
        ccd.u0(ccd.u0() * x_ratio);
        ccd.v0(ccd.v0() * y_ratio);
        ccd.width(new_width);
        ccd.height(new_height);
    }
    
    // access to camera model parameters
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& fx() const { return data()[0]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& fy() const { return data()[1]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& u0() const { return data()[2]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& v0() const { return data()[3]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k1() const { return data()[4]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k2() const { return data()[5]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& p1() const { return data()[6]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& p2() const { return data()[7]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& skew() const { return data()[8]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& width() const { return data()[9]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& height() const { return data()[10]; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fx(const Scalar& v) { data()[0] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fy(const Scalar& v) { data()[1] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void u0(const Scalar& v) { data()[2] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void v0(const Scalar& v) { data()[3] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k1(const Scalar& v) { data()[4] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k2(const Scalar& v) { data()[5] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void p1(const Scalar& v) { data()[6] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void p2(const Scalar& v) { data()[7] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void skew(const Scalar& v) { data()[8] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void width(const Scalar& v) { data()[9] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void height(const Scalar& v) { data()[10] = v; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar* data() { return access_nonconst().data(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar* data() const { return access().data(); }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE PinholeDistanceDistorted<Scalar> getIdeal() const
    {
        return PinholeDistanceDistorted<Scalar>(fx(),fy(),u0(),v0(),width(),height());
    }
    
#ifdef CAMERA_MODELS_HAVE_SERIALIZER
    template <typename Archive>
    void serialize(Archive & archive, std::uint32_t const version)
    { 
        CAMERA_MODELS_SERIALIZE(archive,"fx",data()[0]);
        CAMERA_MODELS_SERIALIZE(archive,"fy",data()[1]);
        CAMERA_MODELS_SERIALIZE(archive,"u0",data()[2]);
        CAMERA_MODELS_SERIALIZE(archive,"v0",data()[3]);
        CAMERA_MODELS_SERIALIZE(archive,"k1",data()[4]);
        CAMERA_MODELS_SERIALIZE(archive,"k2",data()[5]);
        CAMERA_MODELS_SERIALIZE(archive,"p1",data()[6]);
        CAMERA_MODELS_SERIALIZE(archive,"p2",data()[7]);
        CAMERA_MODELS_SERIALIZE(archive,"skew",data()[8]);
        CAMERA_MODELS_SERIALIZE(archive,"width",data()[9]);
        CAMERA_MODELS_SERIALIZE(archive,"height",data()[10]);
    }
#endif // CAMERA_MODELS_HAVE_SERIALIZER

};

/**
 * PinholeDistanceDistorted, Eigen storage.
 */
template<typename _Scalar, int _Options>
class PinholeDistanceDistorted : public PinholeDistanceDistortedBase<PinholeDistanceDistorted<_Scalar, _Options>>
{
    typedef PinholeDistanceDistortedBase<PinholeDistanceDistorted<_Scalar, _Options>> Base;
public:
    typedef typename Eigen::internal::traits<PinholeDistanceDistorted<_Scalar,_Options> >::Scalar Scalar;
    typedef typename Eigen::internal::traits<PinholeDistanceDistorted<_Scalar,_Options> >::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<PinholeDistanceDistorted<_Scalar,_Options> >::ComplexType& ConstComplexReference;
    
    template<class OtherT> using GetOtherType = PinholeDistanceDistorted<OtherT, _Options>;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class PinholeDistanceDistortedBase<PinholeDistanceDistorted<_Scalar, _Options>>;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    EIGEN_DEVICE_FUNC inline PinholeDistanceDistorted(Scalar fx, Scalar fy, Scalar u0, Scalar v0, 
                                                      Scalar k1 = 0.0, Scalar k2 = 0.0, Scalar p1 = 0.0, Scalar p2 = 0.0, 
                                                      Scalar skew = 0.0, Scalar w = 0.0, Scalar h = 0.0)
    {
        access_nonconst() << fx , fy , u0 , v0 , k1 , k2 , p1 , p2 , skew, w, h;
    }
    
    EIGEN_DEVICE_FUNC inline PinholeDistanceDistorted() : parameters(Eigen::Matrix<Scalar,NumParameters,1>::Zero()) { }
    
    EIGEN_DEVICE_FUNC inline PinholeDistanceDistorted(const typename Eigen::internal::traits<PinholeDistanceDistorted<_Scalar,_Options> >::ComplexType& vec) : parameters(vec) 
    {
      
    }
    
    EIGEN_DEVICE_FUNC inline PinholeDistanceDistorted& 
        operator=(const typename Eigen::internal::traits<PinholeDistanceDistorted<_Scalar,_Options> >::ComplexType& vec) 
    { 
        access_nonconst() = vec; 
        return *this; 
    }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    Eigen::Matrix<Scalar,NumParameters,1> parameters;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const PinholeDistanceDistorted<T>& p)
{
    os << "PinholeDistanceDistorted(fx = " << p.fx() << ", fy = " << p.fy() << ", u0 = " << p.u0() << ", v0 = " << p.v0() 
       << ", k1 = " << p.k1() << ", k2 = " << p.k2() << ", p1 = " << p.p1() << ", p2 = " << p.p2() 
       << ", s = " << p.skew() << "," << p.width() << " x " << p.height() << ")";
    return os;
}

}

namespace Eigen 
{
/**
 * PinholeDistanceDistorted, Eigen Map.
 */
template<typename _Scalar, int _Options>
class Map<cammod::PinholeDistanceDistorted<_Scalar>, _Options>
    : public cammod::PinholeDistanceDistortedBase<Map<cammod::PinholeDistanceDistorted<_Scalar>, _Options>>
{
    typedef cammod::PinholeDistanceDistortedBase<Map<cammod::PinholeDistanceDistorted<_Scalar>, _Options>> Base;
    
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef typename internal::traits<Map>::ComplexType& ComplexReference;
    typedef const typename internal::traits<Map>::ComplexType& ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class cammod::PinholeDistanceDistortedBase<Map<cammod::PinholeDistanceDistorted<_Scalar>, _Options>>;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(Scalar* coeffs) : parameters(coeffs)  { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    
    Map<Matrix<Scalar,NumParameters,1>,_Options> parameters;
};

/**
 * PinholeDistanceDistorted, Eigen Map const.
 */
template<typename _Scalar, int _Options>
class Map<const cammod::PinholeDistanceDistorted<_Scalar>, _Options>
    : public cammod::PinholeDistanceDistortedBase<Map<const cammod::PinholeDistanceDistorted<_Scalar>, _Options>>
{
    typedef cammod::PinholeDistanceDistortedBase<Map<const cammod::PinholeDistanceDistorted<_Scalar>, _Options>> Base;
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
    
#endif // PINHOLE_DISTANCE_DISTORTED_CAMERA_MODEL_HPP
