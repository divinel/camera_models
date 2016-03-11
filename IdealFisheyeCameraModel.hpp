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
 * Ideal OpenCV 3.0 Fisheye camera model.
 * ****************************************************************************
 */

#ifndef IDEAL_FISHEYE_CAMERA_MODEL_HPP
#define IDEAL_FISHEYE_CAMERA_MODEL_HPP

#include <CameraModelHelpers.hpp>

// fwd
namespace camera 
{
template<typename _Scalar, int _Options = 0> class IdealFisheyeCameraModel;

namespace internal
{
static constexpr unsigned int IdealFisheyeCameraModelParameterCount = 7;
}
}

// Eigen Traits 
namespace Eigen 
{
    namespace internal 
    {
        template<typename _Scalar, int _Options>
        struct traits<camera::IdealFisheyeCameraModel<_Scalar,_Options> > 
        {
            typedef _Scalar Scalar;
            typedef Matrix<Scalar,camera::internal::IdealFisheyeCameraModelParameterCount,1> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<camera::IdealFisheyeCameraModel<_Scalar>, _Options> > : traits<camera::IdealFisheyeCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<Matrix<Scalar,camera::internal::IdealFisheyeCameraModelParameterCount,1>,_Options> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<const camera::IdealFisheyeCameraModel<_Scalar>, _Options> > : traits<const camera::IdealFisheyeCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<const Matrix<Scalar,camera::internal::IdealFisheyeCameraModelParameterCount,1>,_Options> ComplexType;
        };    
    }
}

namespace camera
{
/**
 * Fisheye Camera Model, Model Specific Functions.
 */
template<typename Derived>
class IdealFisheyeCameraModelBase : public CameraFunctions<Derived>, public ComplexTypes<typename Eigen::internal::traits<Derived>::Scalar>
{
    typedef CameraFunctions<Derived> FunctionsBase;
public:
    // various helpers and types
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
    typedef typename Eigen::internal::traits<Derived>::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Derived>::ComplexType& ConstComplexReference;
    static constexpr CameraModelType ModelType = CameraModelType::IdealFisheye;
    
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
    EIGEN_DEVICE_FUNC inline IdealFisheyeCameraModel<NewScalarType> cast() const { return IdealFisheyeCameraModel<NewScalarType>(access().template cast<NewScalarType>()); }
    
    template<typename OtherDerived> 
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE IdealFisheyeCameraModelBase<Derived>& operator=(const IdealFisheyeCameraModelBase<OtherDerived> & other) { access_nonconst() = other.access(); return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const  { return static_cast<const Derived*>(this)->access(); }
private:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst()  { return static_cast<Derived*>(this)->access_nonconst(); }
public:
        
    static constexpr unsigned int NumParameters = camera::internal::IdealFisheyeCameraModelParameterCount;
    static constexpr unsigned int ParametersToOptimize = NumParameters - 3;
    static constexpr bool CalibrationSupported = true;
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverse(const Derived& ccd, T x, T y) 
    {
        typename ComplexTypes<T>::PointT ret;
        
        typename ComplexTypes<T>::PixelT invKpt( (x - ccd.u0()) / ccd.fx() , (y - ccd.v0()) / ccd.fy() );
        
        T scale = T(1.0);
        
        T theta_d = invKpt.norm();
        if(theta_d > T(1e-8))
        {
            scale = tan(theta_d) / theta_d;
        }
        
        ret(0) = invKpt(0) * scale;
        ret(1) = invKpt(1) * scale;
        ret(2) = T(1.0);
   
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt) 
    {
        typename ComplexTypes<T>::PixelT ret;
        
        const T a = tmp_pt(0) / tmp_pt(2);
        const T b = tmp_pt(1) / tmp_pt(2);
        const T r2 = a*a + b*b;
        const T r = sqrt(r2);
        const T theta = atan(r);
        
        T inv_r = r > T(1e-8) ? T(1.0)/r : T(1.0);
        T cdist = r > T(1e-8) ? theta * inv_r : T(1.0);
        
        const typename ComplexTypes<T>::PixelT xd1(a * cdist, b * cdist);
        const typename ComplexTypes<T>::PixelT xd3(xd1(0), xd1(1));
        ret = typename ComplexTypes<T>::PixelT(xd3(0) * ccd.fx() + ccd.u0(), xd3(1) * ccd.fy() + ccd.v0()); 
        
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
        const T x2 = (x - ccd.u0()) * (x - ccd.u0());
        const T y2 = (y - ccd.v0()) * (y - ccd.v0());
        if(x2 + y2 < ccd.radius() * ccd.radius())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void resizeViewport(Derived& ccd, const T& new_width, const T& new_height)
    {
        const T x_ratio = new_width / ccd.width();
        const T y_ratio = new_height / ccd.height();
        const T r_ratio = std::min(x_ratio, y_ratio);
        
        ccd.fx(ccd.fx() * x_ratio);
        ccd.fy(ccd.fy() * y_ratio);
        ccd.u0(ccd.u0() * x_ratio);
        ccd.v0(ccd.v0() * y_ratio);
        ccd.radius(ccd.radius() * r_ratio);
        ccd.width(new_width);
        ccd.height(new_height);
    }
    
    // access to camera model parameters
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& fx() const { return data()[0]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& fy() const { return data()[1]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& u0() const { return data()[2]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& v0() const { return data()[3]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& width() const { return data()[4]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& height() const { return data()[5]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& radius() const { return data()[6]; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fx(const Scalar& v) { data()[0] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fy(const Scalar& v) { data()[1] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void u0(const Scalar& v) { data()[2] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void v0(const Scalar& v) { data()[3] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void width(const Scalar& v) { data()[4] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void height(const Scalar& v) { data()[5] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void radius(const Scalar& v) { data()[6] = v; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar* data() { return access_nonconst().data(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar* data() const { return access().data(); }
    
#ifdef CAMERA_MODELS_HAVE_SERIALIZER
    template <typename Archive>
    void serialize(Archive & archive, std::uint32_t const version)
    { 
        CAMERA_MODELS_SERIALIZE(archive,"fx",data()[0]);
        CAMERA_MODELS_SERIALIZE(archive,"fy",data()[1]);
        CAMERA_MODELS_SERIALIZE(archive,"u0",data()[2]);
        CAMERA_MODELS_SERIALIZE(archive,"v0",data()[3]);
        CAMERA_MODELS_SERIALIZE(archive,"width",data()[4]);
        CAMERA_MODELS_SERIALIZE(archive,"height",data()[5]);
        CAMERA_MODELS_SERIALIZE(archive,"radius",data()[6]);
    }
#endif // CAMERA_MODELS_HAVE_SERIALIZER
};

/**
 * Fisheye Camera Model, Eigen storage.
 */
template<typename _Scalar, int _Options>
class IdealFisheyeCameraModel : public IdealFisheyeCameraModelBase<IdealFisheyeCameraModel<_Scalar, _Options>>
{
    typedef IdealFisheyeCameraModelBase<IdealFisheyeCameraModel<_Scalar, _Options>> Base;
public:
    typedef typename Eigen::internal::traits<IdealFisheyeCameraModel<_Scalar,_Options> >::Scalar Scalar;
    typedef typename Eigen::internal::traits<IdealFisheyeCameraModel<_Scalar,_Options> >::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<IdealFisheyeCameraModel<_Scalar,_Options> >::ComplexType& ConstComplexReference;
    
    template<class OtherT> using GetOtherType = IdealFisheyeCameraModel<OtherT, _Options>;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class IdealFisheyeCameraModelBase<IdealFisheyeCameraModel<_Scalar, _Options>>;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    EIGEN_DEVICE_FUNC inline IdealFisheyeCameraModel(Scalar fx, Scalar fy, Scalar u0, Scalar v0, Scalar w = 0.0, Scalar h = 0.0, Scalar rad = 0.0)
    {
        access_nonconst() << fx , fy , u0 , v0 , w, h, rad;
    }
    
    EIGEN_DEVICE_FUNC inline IdealFisheyeCameraModel() : parameters(Eigen::Matrix<Scalar,NumParameters,1>::Zero()) { }
    EIGEN_DEVICE_FUNC inline IdealFisheyeCameraModel(const typename Eigen::internal::traits<IdealFisheyeCameraModel<_Scalar,_Options> >::ComplexType& vec) : parameters(vec) { }
    
    EIGEN_DEVICE_FUNC inline IdealFisheyeCameraModel& operator=(const typename Eigen::internal::traits<IdealFisheyeCameraModel<_Scalar,_Options> >::ComplexType& vec) { access_nonconst() = vec; return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    Eigen::Matrix<Scalar,NumParameters,1> parameters;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const IdealFisheyeCameraModel<T>& p)
{
    os << "IdealFisheyeCameraModel(fx = " << p.fx() << ", fy = " << p.fy() << ", u0 = " << p.u0() << ", v0 = " << p.v0() << ", " << p.width() << " x " << p.height() << ", radius = " << p.radius() << ")";
    return os;
}

}

namespace Eigen 
{
/**
 * Fisheye Camera Model, Eigen Map.
 */
template<typename _Scalar, int _Options>
class Map<camera::IdealFisheyeCameraModel<_Scalar>, _Options> : public camera::IdealFisheyeCameraModelBase<Map<camera::IdealFisheyeCameraModel<_Scalar>, _Options>>
{
    typedef camera::IdealFisheyeCameraModelBase<Map<camera::IdealFisheyeCameraModel<_Scalar>, _Options>> Base;
    
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef typename internal::traits<Map>::ComplexType& ComplexReference;
    typedef const typename internal::traits<Map>::ComplexType& ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class camera::IdealFisheyeCameraModelBase<Map<camera::IdealFisheyeCameraModel<_Scalar>, _Options>>;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(Scalar* coeffs) : parameters(coeffs)  { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    
    Map<Matrix<Scalar,NumParameters,1>,_Options> parameters;
};

/**
 * Fisheye Camera Model, Eigen Map const.
 */
template<typename _Scalar, int _Options>
class Map<const camera::IdealFisheyeCameraModel<_Scalar>, _Options> : public camera::IdealFisheyeCameraModelBase<Map<const camera::IdealFisheyeCameraModel<_Scalar>, _Options>>
{
    typedef camera::IdealFisheyeCameraModelBase<Map<const camera::IdealFisheyeCameraModel<_Scalar>, _Options>> Base;
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
    
#endif // IDEAL_FISHEYE_CAMERA_MODEL_HPP
