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
 * Distorted OpenCV 3.0 Fisheye camera model.
 * ****************************************************************************
 */

#ifndef FISHEYE_CAMERA_MODEL_HPP
#define FISHEYE_CAMERA_MODEL_HPP

#include <CameraModelHelpers.hpp>
#include <IdealFisheyeCameraModel.hpp>

// fwd
namespace camera 
{
template<typename _Scalar, int _Options = 0> class FisheyeCameraModel;

namespace internal
{
static constexpr unsigned int FisheyeCameraModelParameterCount = 12;
}
}

// Eigen Traits 
namespace Eigen 
{
    namespace internal 
    {
        template<typename _Scalar, int _Options>
        struct traits<camera::FisheyeCameraModel<_Scalar,_Options> > 
        {
            typedef _Scalar Scalar;
            typedef Matrix<Scalar,camera::internal::FisheyeCameraModelParameterCount,1> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<camera::FisheyeCameraModel<_Scalar>, _Options> > : traits<camera::FisheyeCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<Matrix<Scalar,camera::internal::FisheyeCameraModelParameterCount,1>,_Options> ComplexType;
        };
        
        template<typename _Scalar, int _Options>
        struct traits<Map<const camera::FisheyeCameraModel<_Scalar>, _Options> > : traits<const camera::FisheyeCameraModel<_Scalar, _Options> > 
        {
            typedef _Scalar Scalar;
            typedef Map<const Matrix<Scalar,camera::internal::FisheyeCameraModelParameterCount,1>,_Options> ComplexType;
        };    
    }
}

namespace camera
{
/**
 * Fisheye Camera Model, Model Specific Functions.
 */
template<typename Derived>
class FisheyeCameraModelBase : public CameraFunctions<Derived>, public ComplexTypes<typename Eigen::internal::traits<Derived>::Scalar>
{
    typedef CameraFunctions<Derived> FunctionsBase;
public:
    // various helpers and types
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
    typedef typename Eigen::internal::traits<Derived>::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<Derived>::ComplexType& ConstComplexReference;
    static constexpr CameraModelType ModelType = CameraModelType::Fisheye;
    
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
    EIGEN_DEVICE_FUNC inline FisheyeCameraModel<NewScalarType> cast() const { return FisheyeCameraModel<NewScalarType>(access().template cast<NewScalarType>()); }
    
    template<typename OtherDerived> 
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE FisheyeCameraModelBase<Derived>& operator=(const FisheyeCameraModelBase<OtherDerived> & other) { access_nonconst() = other.access(); return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const  { return static_cast<const Derived*>(this)->access(); }
private:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst()  { return static_cast<Derived*>(this)->access_nonconst(); }
public:
        
    static constexpr unsigned int NumParameters = camera::internal::FisheyeCameraModelParameterCount;
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
            T theta = theta_d;
            for(unsigned int j = 0; j < 10; ++j)
            {
                T theta2 = theta*theta, 
                  theta4 = theta2*theta2, 
                  theta6 = theta4*theta2, 
                  theta8 = theta6*theta2;
                theta = theta_d / (T(1.0) + ccd.k1() * theta2 + ccd.k2() * theta4 + ccd.k3() * theta6 + ccd.k4() * theta8);
            }
            
            scale = tan(theta) / theta_d;
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
        
        const T theta2 = theta*theta, 
                theta3 = theta2*theta, 
                theta4 = theta2*theta2, 
                theta5 = theta4*theta,
                theta6 = theta3*theta3, 
                theta7 = theta6*theta, 
                theta8 = theta4*theta4, 
                theta9 = theta8*theta;
        
        const T theta_d = theta + ccd.k1() * theta3 + ccd.k2() * theta5 + ccd.k3()*theta7 + ccd.k4() * theta9;
        
        T inv_r = r > T(1e-8) ? T(1.0)/r : T(1.0);
        T cdist = r > T(1e-8) ? theta_d * inv_r : T(1.0);
        
        const typename ComplexTypes<T>::PixelT xd1(a * cdist, b * cdist);
        const typename ComplexTypes<T>::PixelT xd3(xd1(0) + ccd.skew() * xd1(1), xd1(1));
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
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k1() const { return data()[4]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k2() const { return data()[5]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k3() const { return data()[6]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& k4() const { return data()[7]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& skew() const { return data()[8]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& width() const { return data()[9]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& height() const { return data()[10]; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar& radius() const { return data()[11]; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fx(const Scalar& v) { data()[0] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void fy(const Scalar& v) { data()[1] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void u0(const Scalar& v) { data()[2] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void v0(const Scalar& v) { data()[3] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k1(const Scalar& v) { data()[4] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k2(const Scalar& v) { data()[5] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k3(const Scalar& v) { data()[6] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void k4(const Scalar& v) { data()[7] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void skew(const Scalar& v) { data()[8] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void width(const Scalar& v) { data()[9] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void height(const Scalar& v) { data()[10] = v; }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void radius(const Scalar& v) { data()[11] = v; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar* data() { return access_nonconst().data(); }
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar* data() const { return access().data(); }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE IdealFisheyeCameraModel<Scalar> getIdeal() const
    {
        return IdealFisheyeCameraModel<Scalar>(fx(),fy(),u0(),v0(),width(),height(),radius());
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
        CAMERA_MODELS_SERIALIZE(archive,"k3",data()[6]);
        CAMERA_MODELS_SERIALIZE(archive,"k4",data()[7]);
        CAMERA_MODELS_SERIALIZE(archive,"skew",data()[8]);
        CAMERA_MODELS_SERIALIZE(archive,"width",data()[9]);
        CAMERA_MODELS_SERIALIZE(archive,"height",data()[10]);
        CAMERA_MODELS_SERIALIZE(archive,"radius",data()[11]);
    }
#endif // CAMERA_MODELS_HAVE_SERIALIZER
};

/**
 * Fisheye Camera Model, Eigen storage.
 */
template<typename _Scalar, int _Options>
class FisheyeCameraModel : public FisheyeCameraModelBase<FisheyeCameraModel<_Scalar, _Options>>
{
    typedef FisheyeCameraModelBase<FisheyeCameraModel<_Scalar, _Options>> Base;
public:
    typedef typename Eigen::internal::traits<FisheyeCameraModel<_Scalar,_Options> >::Scalar Scalar;
    typedef typename Eigen::internal::traits<FisheyeCameraModel<_Scalar,_Options> >::ComplexType& ComplexReference;
    typedef const typename Eigen::internal::traits<FisheyeCameraModel<_Scalar,_Options> >::ComplexType& ConstComplexReference;
    
    template<class OtherT> using GetOtherType = FisheyeCameraModel<OtherT, _Options>;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class FisheyeCameraModelBase<FisheyeCameraModel<_Scalar, _Options>>;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    EIGEN_DEVICE_FUNC inline FisheyeCameraModel(Scalar fx, Scalar fy, Scalar u0, Scalar v0, Scalar k1 = 0.0, Scalar k2 = 0.0, Scalar k3 = 0.0, Scalar k4 = 0.0, Scalar skew = 0.0, Scalar w = 0.0, Scalar h = 0.0, Scalar rad = 0.0)
    {
        access_nonconst() << fx , fy , u0 , v0 , k1 , k2 , k3 , k4 , skew , w , h, rad;
    }
    
    EIGEN_DEVICE_FUNC inline FisheyeCameraModel() : parameters(Eigen::Matrix<Scalar,NumParameters,1>::Zero()) { }
    EIGEN_DEVICE_FUNC inline FisheyeCameraModel(const typename Eigen::internal::traits<FisheyeCameraModel<_Scalar,_Options> >::ComplexType& vec) : parameters(vec) { }
    
    EIGEN_DEVICE_FUNC inline FisheyeCameraModel& operator=(const typename Eigen::internal::traits<FisheyeCameraModel<_Scalar,_Options> >::ComplexType& vec) { access_nonconst() = vec; return *this; }
    
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ConstComplexReference access() const { return parameters; }
protected:
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE ComplexReference access_nonconst() { return parameters; }
    Eigen::Matrix<Scalar,NumParameters,1> parameters;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const FisheyeCameraModel<T>& p)
{
    os << "FisheyeCameraModel(fx = " << p.fx() << ", fy = " << p.fy() << ", u0 = " << p.u0() << ", v0 = " << p.v0() << ", k1 = " << p.k1() << ", k2 = " << p.k2() << ", k3 = " << p.k3() << ", k4 = " << p.k4() << ", s = " << p.skew() << ", " << p.width() << " x " << p.height() << ", radius = " << p.radius() << ")";
    return os;
}

}

namespace Eigen 
{
/**
 * Fisheye Camera Model, Eigen Map.
 */
template<typename _Scalar, int _Options>
class Map<camera::FisheyeCameraModel<_Scalar>, _Options> : public camera::FisheyeCameraModelBase<Map<camera::FisheyeCameraModel<_Scalar>, _Options>>
{
    typedef camera::FisheyeCameraModelBase<Map<camera::FisheyeCameraModel<_Scalar>, _Options>> Base;
    
public:
    typedef typename internal::traits<Map>::Scalar Scalar;
    typedef typename internal::traits<Map>::ComplexType& ComplexReference;
    typedef const typename internal::traits<Map>::ComplexType& ConstComplexReference;
    
    static constexpr unsigned int NumParameters = Base::NumParameters;
    static constexpr unsigned int ParametersToOptimize = Base::ParametersToOptimize;
    static constexpr bool CalibrationSupported = Base::CalibrationSupported;
    
    friend class camera::FisheyeCameraModelBase<Map<camera::FisheyeCameraModel<_Scalar>, _Options>>;
    
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
class Map<const camera::FisheyeCameraModel<_Scalar>, _Options> : public camera::FisheyeCameraModelBase<Map<const camera::FisheyeCameraModel<_Scalar>, _Options>>
{
    typedef camera::FisheyeCameraModelBase<Map<const camera::FisheyeCameraModel<_Scalar>, _Options>> Base;
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
    
#endif // FISHEYE_CAMERA_MODEL_HPP
