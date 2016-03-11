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
 * Neither the name of [project] nor the names of its
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
 * Ideal Pinhole Camera model with inverse in image planes.
 * ****************************************************************************
 */

#ifndef PINHOLE_DISPARITY_CAMERA_MODEL_HPP
#define PINHOLE_DISPARITY_CAMERA_MODEL_HPP

#include <CameraModelHelpers.hpp>

// fwd
namespace camera 
{
template<typename _Scalar, int _Options = 0> class PinholeDisparityCameraModel;

namespace internal
{
static constexpr unsigned int PinholeDisparityModelParameterCount = 6;
}
}

// Eigen Traits 
namespace Eigen 
{
    namespace internal 
    {
        template<typename _Scalar, int _Options>
        struct traits<camera::PinholeDisparityCameraModel<_Scalar,_Options> > 
        {
            typedef _Scalar Scalar;
            typedef Matrix<Scalar,camera::internal::PinholeDisparityModelParameterCount,1> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<camera::PinholeDisparityCameraModel<_Scalar>, _Options> > : traits<camera::PinholeDisparityCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<Matrix<Scalar,camera::internal::PinholeDisparityModelParameterCount,1>,_Options> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<const camera::PinholeDisparityCameraModel<_Scalar>, _Options> > : traits<const camera::PinholeDisparityCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<const Matrix<Scalar,camera::internal::PinholeDisparityModelParameterCount,1>,_Options> ComplexType;
        };    
    }
}

namespace camera
{
/**
 * Pinhole Disparity Camera Model, Model Specific Functions.
 */
template<typename Derived>
class PinholeDisparityCameraModelBase : public CameraFunctions<Derived>, public ComplexTypes<typename Eigen::internal::traits<Derived>::Scalar>
{
    typedef CameraFunctions<Derived> FunctionsBase;
public:
    // various helpers and types
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
    typedef typename Eigen::internal::traits<Derived>::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Derived>::ComplexType& ConstComplexReference;
    static constexpr CameraModelType ModelType = CameraModelType::PinholeDisparity;
    
    using FunctionsBase::forward;
    using FunctionsBase::inverse;
    using FunctionsBase::twoFrameProject;
    using FunctionsBase::worldToCamera;
    using FunctionsBase::cameraToWorld;
    using FunctionsBase::pixelValid;
    using FunctionsBase::pixelValidSquare;
    using FunctionsBase::pixelValidCircular;
    using FunctionsBase::resizeViewport;
    
    template<typename NewScalarType>
    EIGEN_DEVICE_FUNC inline PinholeDisparityCameraModel<NewScalarType> cast() const { return PinholeDisparityCameraModel<NewScalarType>(access().template cast<NewScalarType>()); }
    
    template<typename OtherDerived> 
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE PinholeDisparityCameraModelBase<Derived>& operator=(const PinholeDisparityCameraModelBase<OtherDerived> & other) { access_nonconst() = other.access(); return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const  { return static_cast<const Derived*>(this)->access(); }
private:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst()  { return static_cast<Derived*>(this)->access_nonconst(); }
public:
        
    static constexpr unsigned int NumParameters = camera::internal::PinholeDisparityModelParameterCount;
    static constexpr unsigned int ParametersToOptimize = NumParameters - 2;
    static constexpr bool CalibrationSupported = true;
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverse(const Derived& ccd, T x, T y) 
    {
        typename ComplexTypes<T>::PointT ret;

        ret(0) = (x - ccd.u0())/ccd.fx();
        ret(1) = (y - ccd.v0())/ccd.fy();
        ret(2) = T(1.0);
        
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const Derived& ccd, const typename ComplexTypes<T>::PointT& tmp_pt) 
    {
        typename ComplexTypes<T>::PixelT ret;

        ret(0) = (ccd.fx() * tmp_pt(0))/tmp_pt(2) + ccd.u0();
        ret(1) = (ccd.fy() * tmp_pt(1))/tmp_pt(2) + ccd.v0();
        
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
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar getFieldOfViewX() const { return getFieldOfView(fx(),width()); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar getFieldOfViewY() const { return getFieldOfView(fy(),height()); }

    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void setFromFieldOfView(const Scalar& fovx, const Scalar& fovy)
    {
        fx(getFocalLength(fovx, width()));
        fy(getFocalLength(fovy, height()));
    }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Eigen::Matrix<Scalar,3,3> getIntrinsicMatrix() const
    {
        Eigen::Matrix<Scalar,3,3> ret;
        ret << fx() , Scalar(0.0) , u0(),
               Scalar(0.0) , fy() , v0(),
               Scalar(0.0), Scalar(0.0), Scalar(1.0);
        return ret;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void resizeViewport(Derived& ccd, const T& new_width, const T& new_height)
    {
        using std::min;
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
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& width() const { return data()[4]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& height() const { return data()[5]; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fx(const Scalar& v) { data()[0] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fy(const Scalar& v) { data()[1] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void u0(const Scalar& v) { data()[2] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void v0(const Scalar& v) { data()[3] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void width(const Scalar& v) { data()[4] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void height(const Scalar& v) { data()[5] = v; }
    
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
    }
#endif // CAMERA_MODELS_HAVE_SERIALIZER
};

/**
 * Pinhole Disparity Camera Model, Eigen storage.
 */
template<typename _Scalar, int _Options>
class PinholeDisparityCameraModel : public PinholeDisparityCameraModelBase<PinholeDisparityCameraModel<_Scalar, _Options>>
{
    typedef PinholeDisparityCameraModelBase<PinholeDisparityCameraModel<_Scalar, _Options>> Base;
public:
    typedef typename Eigen::internal::traits<PinholeDisparityCameraModel<_Scalar,_Options> >::Scalar Scalar;
    typedef typename Eigen::internal::traits<PinholeDisparityCameraModel<_Scalar,_Options> >::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<PinholeDisparityCameraModel<_Scalar,_Options> >::ComplexType& ConstComplexReference;
    
    template<class OtherT> using GetOtherType = PinholeDisparityCameraModel<OtherT, _Options>;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class PinholeDisparityCameraModelBase<PinholeDisparityCameraModel<_Scalar, _Options>>;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    EIGEN_DEVICE_FUNC inline PinholeDisparityCameraModel(Scalar fx, Scalar fy, Scalar u0, Scalar v0, Scalar w = 0.0, Scalar h = 0.0)
    {
        parameters << fx , fy , u0 , v0 , w , h;
    }
    EIGEN_DEVICE_FUNC inline PinholeDisparityCameraModel() : parameters(Eigen::Matrix<Scalar,NumParameters,1>::Zero()) { }
    EIGEN_DEVICE_FUNC inline PinholeDisparityCameraModel(const typename Eigen::internal::traits<PinholeDisparityCameraModel<_Scalar,_Options> >::ComplexType& vec) : parameters(vec) { }
    
    EIGEN_DEVICE_FUNC inline PinholeDisparityCameraModel& operator=(const typename Eigen::internal::traits<PinholeDisparityCameraModel<_Scalar,_Options> >::ComplexType& vec) { access_nonconst() = vec; return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    Eigen::Matrix<Scalar,NumParameters,1> parameters;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const PinholeDisparityCameraModel<T>& p)
{
    os << "PinholeDisparityCamera(fx = " << p.fx() << ", fy = " << p.fy() << ", u0 = " << p.u0() << ", v0 = " << p.v0() << "," << p.width() << " x " << p.height() << ")";
    return os;
}

}

namespace Eigen 
{
/**
 * Pinhole Disparity Camera Model, Eigen Map.
 */
template<typename _Scalar, int _Options>
class Map<camera::PinholeDisparityCameraModel<_Scalar>, _Options> : public camera::PinholeDisparityCameraModelBase<Map<camera::PinholeDisparityCameraModel<_Scalar>, _Options>>
{
    typedef camera::PinholeDisparityCameraModelBase<Map<camera::PinholeDisparityCameraModel<_Scalar>, _Options>> Base;
    
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef typename internal::traits<Map>::ComplexType& ComplexReference;
    typedef const typename internal::traits<Map>::ComplexType& ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class camera::PinholeDisparityCameraModelBase<Map<camera::PinholeDisparityCameraModel<_Scalar>, _Options>>;
    
    EIGEN_INHERIT_ASSIGNMENT_EQUAL_OPERATOR(Map)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Map(Scalar* coeffs) : parameters(coeffs)  { }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    
    Map<Matrix<Scalar,NumParameters,1>,_Options> parameters;
};

/**
 * Pinhole Disparity Camera Model, Eigen Map const.
 */
template<typename _Scalar, int _Options>
class Map<const camera::PinholeDisparityCameraModel<_Scalar>, _Options> : public camera::PinholeDisparityCameraModelBase<Map<const camera::PinholeDisparityCameraModel<_Scalar>, _Options>>
{
    typedef camera::PinholeDisparityCameraModelBase<Map<const camera::PinholeDisparityCameraModel<_Scalar>, _Options>> Base;
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
    
    #endif // PINHOLE_DISPARITY_CAMERA_MODEL_HPP
