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
 * Spherical Panorama camera model as in PovRay.
 * ****************************************************************************
 */

#ifndef SPHERICAL_POVRAY_CAMERA_MODEL_HPP
#define SPHERICAL_POVRAY_CAMERA_MODEL_HPP

#include <CameraModelHelpers.hpp>

// fwd
namespace camera 
{
template<typename _Scalar, int _Options = 0> class SphericalPovRayCameraModel;

namespace internal
{
static constexpr unsigned int SphericalPovRayModelParameterCount = 4;
}
}

// Eigen Traits 
namespace Eigen 
{
    namespace internal 
    {
        template<typename _Scalar, int _Options>
        struct traits<camera::SphericalPovRayCameraModel<_Scalar,_Options> > 
        {
            typedef _Scalar Scalar;
            typedef Matrix<Scalar,camera::internal::SphericalPovRayModelParameterCount,1> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<camera::SphericalPovRayCameraModel<_Scalar>, _Options> > : traits<camera::SphericalPovRayCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<Matrix<Scalar,camera::internal::SphericalPovRayModelParameterCount,1>,_Options> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<const camera::SphericalPovRayCameraModel<_Scalar>, _Options> > : traits<const camera::SphericalPovRayCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<const Matrix<Scalar,camera::internal::SphericalPovRayModelParameterCount,1>,_Options> ComplexType;
        };    
    }
}

namespace camera
{
/**
 * Spherical PovRay Camera Model, Model Specific Functions.
 */
template<typename Derived>
class SphericalPovRayCameraModelBase : public CameraFunctions<Derived>, public ComplexTypes<typename Eigen::internal::traits<Derived>::Scalar>
{
    typedef CameraFunctions<Derived> FunctionsBase;
public:
    // various helpers and types
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
    typedef typename Eigen::internal::traits<Derived>::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Derived>::ComplexType& ConstComplexReference;
    static constexpr CameraModelType ModelType = CameraModelType::SphericalPovRay;
    
    using FunctionsBase::forward;
    using FunctionsBase::inverse;
    using FunctionsBase::inverseAtDistance;
    using FunctionsBase::twoFrameProject;
    using FunctionsBase::worldToCamera;
    using FunctionsBase::cameraToWorld;
    using FunctionsBase::pixelValid;
    using FunctionsBase::pixelValidSquare;
    using FunctionsBase::pixelValidCircular;
    using FunctionsBase::resizeViewport;
    
    template<typename NewScalarType>
    EIGEN_DEVICE_FUNC inline SphericalPovRayCameraModel<NewScalarType> cast() const { return SphericalPovRayCameraModel<NewScalarType>(access().template cast<NewScalarType>()); }
    
    template<typename OtherDerived> 
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE SphericalPovRayCameraModelBase<Derived>& operator=(const SphericalPovRayCameraModelBase<OtherDerived> & other) { access_nonconst() = other.access(); return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const  { return static_cast<const Derived*>(this)->access(); }
private:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst()  { return static_cast<Derived*>(this)->access_nonconst(); }
public:
        
    static constexpr unsigned int NumParameters = camera::internal::SphericalPovRayModelParameterCount;
    static constexpr bool CalibrationSupported = false;
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverse(const Derived& ccd, T x, T y) 
    {
        typename ComplexTypes<T>::PointT ret;

        const T angleV = ccd.min_angle() + ( (ccd.max_angle() - ccd.min_angle()) / ccd.height () * (y + T(1.0f)));
        const T angleH = (T(M_PI) * ccd.width() - T(2.0f * M_PI) * x) / ccd.width();
        
        // odd PovRay spherical coordinates
        ret(0) = -cos(angleV) * sin(angleH);
        ret(1) = -sin(angleV);
        ret(2) =  cos(angleV) * cos(angleH);

        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt) 
    {
        typename ComplexTypes<T>::PixelT ret;

        const T radius = sqrt(tmp_pt(0) * tmp_pt(0) + tmp_pt(1) * tmp_pt(1) + tmp_pt(2) * tmp_pt(2));
        const T angle1 = -asin(tmp_pt(1) / radius);
        const T angle2 = (atan2(tmp_pt(0), tmp_pt(2)) + T(M_PI));
        
        ret(0) = ccd.width() * ((angle2)/(T(2.0f * M_PI)));
        ret(1) = (ccd.height() * angle1 - (ccd.height() * ccd.min_angle() - (ccd.max_angle() - ccd.min_angle()) )) / (ccd.max_angle() - ccd.min_angle());
        
        if(ret(0) > ccd.width() - T(1.0f))
        {
            ret(0) = T(0.0f);
        }
        
        return ret;
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
 * Spherical PovRay Camera Model, Eigen storage.
 */
template<typename _Scalar, int _Options>
class SphericalPovRayCameraModel : public SphericalPovRayCameraModelBase<SphericalPovRayCameraModel<_Scalar, _Options>>
{
    typedef SphericalPovRayCameraModelBase<SphericalPovRayCameraModel<_Scalar, _Options>> Base;
public:
    typedef typename Eigen::internal::traits<SphericalPovRayCameraModel<_Scalar,_Options> >::Scalar Scalar;
    typedef typename Eigen::internal::traits<SphericalPovRayCameraModel<_Scalar,_Options> >::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<SphericalPovRayCameraModel<_Scalar,_Options> >::ComplexType& ConstComplexReference;
    
    template<class OtherT> using GetOtherType = SphericalPovRayCameraModel<OtherT, _Options>;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class SphericalPovRayCameraModelBase<SphericalPovRayCameraModel<_Scalar, _Options>>;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    EIGEN_DEVICE_FUNC inline SphericalPovRayCameraModel(Scalar w, Scalar h, Scalar min_a, Scalar max_a)
    {
        access_nonconst() << w , h , min_a , max_a;
    }
    EIGEN_DEVICE_FUNC inline SphericalPovRayCameraModel() : parameters(Eigen::Matrix<Scalar,NumParameters,1>::Zero()) { }
    EIGEN_DEVICE_FUNC inline SphericalPovRayCameraModel(const typename Eigen::internal::traits<SphericalPovRayCameraModel<_Scalar,_Options> >::ComplexType& vec) : parameters(vec) { }
    
    EIGEN_DEVICE_FUNC inline SphericalPovRayCameraModel& operator=(const typename Eigen::internal::traits<SphericalPovRayCameraModel<_Scalar,_Options> >::ComplexType& vec) { access_nonconst() = vec; return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    Eigen::Matrix<Scalar,NumParameters,1> parameters;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const SphericalPovRayCameraModel<T>& p)
{
    os << "SphericalPovRayCameraModel(width = " << p.width() << ", height = " << p.height() << ", min_angle = " << p.min_angle() << ", max_angle = " << p.max_angle() << ")";
    return os;
}

}

namespace Eigen 
{
/**
 * Spherical PovRay Camera Model, Eigen Map.
 */
template<typename _Scalar, int _Options>
class Map<camera::SphericalPovRayCameraModel<_Scalar>, _Options> : public camera::SphericalPovRayCameraModelBase<Map<camera::SphericalPovRayCameraModel<_Scalar>, _Options>>
{
    typedef camera::SphericalPovRayCameraModelBase<Map<camera::SphericalPovRayCameraModel<_Scalar>, _Options>> Base;
    
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef typename internal::traits<Map>::ComplexType& ComplexReference;
    typedef const typename internal::traits<Map>::ComplexType& ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class camera::SphericalPovRayCameraModelBase<Map<camera::SphericalPovRayCameraModel<_Scalar>, _Options>>;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(Scalar* coeffs) : parameters(coeffs)  { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    
    Map<Matrix<Scalar,NumParameters,1>,_Options> parameters;
};

/**
 * Spherical PovRay Camera Model, Eigen Map const.
 */
template<typename _Scalar, int _Options>
class Map<const camera::SphericalPovRayCameraModel<_Scalar>, _Options> : public camera::SphericalPovRayCameraModelBase<Map<const camera::SphericalPovRayCameraModel<_Scalar>, _Options>>
{
    typedef camera::SphericalPovRayCameraModelBase<Map<const camera::SphericalPovRayCameraModel<_Scalar>, _Options>> Base;
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
    
#endif // SPHERICAL_POVRAY_CAMERA_MODEL_HPP
