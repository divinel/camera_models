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
 * Simple Spherical Panorama camera model.
 * ****************************************************************************
 */

#ifndef SPHERICAL_CAMERA_MODEL_HPP
#define SPHERICAL_CAMERA_MODEL_HPP

#include <CameraModels/CameraModelUtils.hpp>

// fwd
namespace cammod 
{
template<typename _Scalar, int _Options = 0> class Spherical;

namespace internal
{
static constexpr unsigned int SphericalModelParameterCount = 4;
}
}

// Eigen Traits, but also some traits for our use
namespace Eigen 
{
    namespace internal 
    {
        template<typename _Scalar, int _Options>
        struct traits<cammod::Spherical<_Scalar,_Options> > 
        {
            typedef _Scalar Scalar;
            typedef Matrix<Scalar,cammod::internal::SphericalModelParameterCount,1> ComplexType;
            static constexpr bool HasForwardPointJacobian = false;
            static constexpr bool HasForwardParametersJacobian = false;
            static constexpr bool HasInversePointJacobian = false;
            static constexpr bool HasInverseParametersJacobian = false;
            static constexpr unsigned int NumParameters = cammod::internal::SphericalModelParameterCount;
            static constexpr unsigned int ParametersToOptimize = 0;
            static constexpr bool CalibrationSupported = false;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<cammod::Spherical<_Scalar>, _Options> >
            : traits<cammod::Spherical<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<Matrix<Scalar,cammod::internal::SphericalModelParameterCount,1>,_Options> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<const cammod::Spherical<_Scalar>, _Options> >
            : traits<const cammod::Spherical<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<const Matrix<Scalar,cammod::internal::SphericalModelParameterCount,1>,_Options> ComplexType;
        };    
    }
}

namespace cammod
{
/**
 * Spherical Camera Model, Model Specific Functions.
 */
template<typename Derived>
class SphericalBase : public CameraFunctions<Derived>, 
                      public ComplexTypes<typename Eigen::internal::traits<Derived>::Scalar>
{
    typedef CameraFunctions<Derived> FunctionsBase;
public:
    // various helpers and types
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
    typedef typename Eigen::internal::traits<Derived>::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Derived>::ComplexType& ConstComplexReference;
    static constexpr CameraModelType ModelType = CameraModelType::Spherical;
    
    using FunctionsBase::forward;
    using FunctionsBase::inverse;
    using FunctionsBase::inverseAtDistance;
    using FunctionsBase::twoFrameProject;
    using FunctionsBase::worldToCamera;
    using FunctionsBase::cameraToWorld;
    using FunctionsBase::pointValid;
    using FunctionsBase::pixelValid;
    using FunctionsBase::pixelValidSquare;
    using FunctionsBase::pixelValidCircular;
    using FunctionsBase::resizeViewport;
    
    template<typename NewScalarType>
    EIGEN_DEVICE_FUNC inline Spherical<NewScalarType> cast() const 
    { 
        return Spherical<NewScalarType>(access().template cast<NewScalarType>()); 
    }
    
    template<typename OtherDerived> 
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE SphericalBase<Derived>& 
        operator=(const SphericalBase<OtherDerived> & other) 
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
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT 
        inverse(const Derived& ccd, T x, T y) 
    {
        using Eigen::numext::sin;
        using Eigen::numext::cos;
        
        typename ComplexTypes<T>::PointT ret;

        // camera coordinate frame
        const T angle_horiz = T(M_PI) - ((x / ccd.width()) * T(2.0f * M_PI));
        const T angle_vert = ccd.min_angle() + ((y / ccd.height()) * (ccd.max_angle() - ccd.min_angle()));
        
        // this is a unit vector?
        ret(0) = T(1.0f) * sin(angle_vert) * cos(angle_horiz);
        ret(1) = T(1.0f) * sin(angle_vert) * sin(angle_horiz);
        ret(2) = T(1.0f) * cos(angle_vert);

        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT 
        forward(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt) 
    {
        using Eigen::numext::sqrt;
        using Eigen::numext::acos;
        using Eigen::numext::atan2;
        
        typename ComplexTypes<T>::PixelT ret;

        // camera coordinate frame - convert to spherical coordinates
        const T radius = sqrt(tmp_pt(0) * tmp_pt(0) + tmp_pt(1) * tmp_pt(1) + tmp_pt(2) * tmp_pt(2));
        const T angle1 = acos(tmp_pt(2) / radius);
        const T angle2 = atan2(tmp_pt(1), tmp_pt(0)) + T(M_PI);
        
        // and now to pixels
        ret(0) = ccd.width() - (angle2 / T(2.0f * M_PI)) * ccd.width(); // x
        ret(1) = ((angle1 - ccd.min_angle()) / (ccd.max_angle() - ccd.min_angle())) * ccd.height(); // y
        
        // return
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pointValid(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt)
    {
        return true;
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
        ccd.width(new_width);
        ccd.height(new_height);
    }
    
    // access to camera model parameters
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& width() const { return data()[0]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& height() const { return data()[1]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& min_angle() const { return data()[2]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& max_angle() const { return data()[3]; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void width(const Scalar& v) { data()[0] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void height(const Scalar& v) { data()[1] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void min_angle(const Scalar& v) { data()[2] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void max_angle(const Scalar& v) { data()[3] = v; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar* data() { return access_nonconst().data(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar* data() const { return access().data(); }
    
#ifdef CAMERA_MODELS_HAVE_SERIALIZER
    template <typename Archive>
    void serialize(Archive & archive, std::uint32_t const version)
    { 
        CAMERA_MODELS_SERIALIZE(archive,"width",data()[0]);
        CAMERA_MODELS_SERIALIZE(archive,"height",data()[1]);
        CAMERA_MODELS_SERIALIZE(archive,"min_angle",data()[2]);
        CAMERA_MODELS_SERIALIZE(archive,"max_angle",data()[3]);
    }
#endif // CAMERA_MODELS_HAVE_SERIALIZER
};

/**
 * Spherical Camera Model, Eigen storage.
 */
template<typename _Scalar, int _Options>
class Spherical : public SphericalBase<Spherical<_Scalar, _Options>>
{
    typedef SphericalBase<Spherical<_Scalar, _Options>> Base;
public:
    typedef typename Eigen::internal::traits<Spherical<_Scalar,_Options> >::Scalar Scalar;
    typedef typename Eigen::internal::traits<Spherical<_Scalar,_Options> >::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Spherical<_Scalar,_Options> >::ComplexType& ConstComplexReference;
    
    template<class OtherT> using GetOtherType = Spherical<OtherT, _Options>;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class SphericalBase<Spherical<_Scalar, _Options>>;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    EIGEN_DEVICE_FUNC inline Spherical(Scalar w, Scalar h, Scalar min_a, Scalar max_a)
    {
        access_nonconst() << w , h , min_a , max_a;
    }
    
    inline Spherical() : parameters(Eigen::Matrix<Scalar,NumParameters,1>::Zero()) { }
    
    EIGEN_DEVICE_FUNC inline Spherical(const typename Eigen::internal::traits<Spherical<_Scalar,_Options> >::ComplexType& vec)
      : parameters(vec) 
    { 
      
    }
    
    EIGEN_DEVICE_FUNC inline Spherical& 
        operator=(const typename Eigen::internal::traits<Spherical<_Scalar,_Options> >::ComplexType& vec) 
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
inline std::ostream& operator<<(std::ostream& os, const Spherical<T>& p)
{
    os << "Spherical(width = " << p.width() << ", height = " << p.height() 
       << ", min_angle = " << p.min_angle() << ", max_angle = " << p.max_angle() << ")";
    return os;
}

}

namespace Eigen 
{
/**
  * Spherical Camera Model, Eigen Map.
  */
template<typename _Scalar, int _Options>
class Map<cammod::Spherical<_Scalar>, _Options>
    : public cammod::SphericalBase<Map<cammod::Spherical<_Scalar>, _Options>>
{
    typedef cammod::SphericalBase<Map<cammod::Spherical<_Scalar>, _Options>> Base;
    
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef typename internal::traits<Map>::ComplexType& ComplexReference;
    typedef const typename internal::traits<Map>::ComplexType& ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class cammod::SphericalBase<Map<cammod::Spherical<_Scalar>, _Options>>;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(Scalar* coeffs) : parameters(coeffs)  { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    
    Map<Matrix<Scalar,NumParameters,1>,_Options> parameters;
};

/**
  * Spherical Camera Model, Eigen Map const.
  */
template<typename _Scalar, int _Options>
class Map<const cammod::Spherical<_Scalar>, _Options>
    : public cammod::SphericalBase<Map<const cammod::Spherical<_Scalar>, _Options>>
{
    typedef cammod::SphericalBase<Map<const cammod::Spherical<_Scalar>, _Options>> Base;
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef const typename internal::traits<Map>::ComplexType & ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(const Scalar* coeffs) : parameters(coeffs) { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const  { return parameters; }
protected:
    const Map<const Matrix<Scalar,NumParameters,1>,_Options> parameters;
};

}
    
#endif // SPHERICAL_CAMERA_MODEL_HPP
