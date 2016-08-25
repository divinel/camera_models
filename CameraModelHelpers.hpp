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
 * Common types and functions.
 * ****************************************************************************
 */

#ifndef CAMERA_MODEL_HELPERS_HPP
#define CAMERA_MODEL_HELPERS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <sophus/so2.hpp>
#include <sophus/se2.hpp>
#include <sophus/so3.hpp>
#include <sophus/se3.hpp>

// For future Eigen + CUDA
#ifndef EIGEN_DEVICE_FUNC
#define EIGEN_DEVICE_FUNC
#endif // EIGEN_DEVICE_FUNC

// If Cereal serializer is preferred
#ifdef CAMERA_MODELS_SERIALIZER_CEREAL
#define CAMERA_MODELS_HAVE_SERIALIZER
#define CAMERA_MODELS_SERIALIZE(ARCHIVE,NAME,VAR) ARCHIVE(cereal::make_nvp(NAME,VAR))
#endif // CAMERA_MODELS_SERIALIZER_CEREAL

// If Boost serializer is preferred
#ifdef CAMERA_MODELS_SERIALIZER_BOOST
#define CAMERA_MODELS_HAVE_SERIALIZER
#define CAMERA_MODELS_SERIALIZE(ARCHIVE,NAME,VAR) ARCHIVE & boost::serialization::make_nvp(NAME,VAR)
#endif // CAMERA_MODELS_SERIALIZER_BOOST 

namespace camera
{

/**
 * Camera Models
 */
enum class CameraModelType
{
    Pinhole = 0,
    PinholeDistorted,
    IdealGeneric,
    FullGeneric,
    Spherical,
    SphericalPovRay,
    Fisheye,
    IdealFisheye,
    PinholeDisparity,
    PinholeDisparityDistorted,
    PinholeDisparityBrownConrady
};

template<CameraModelType cmt>
struct CameraModelToTypeAndName;

/**
 * Collection of common 2D/3D types.
 */
template<typename T>
struct ComplexTypes
{
    typedef Sophus::SO3Group<T> RotationT;
    typedef Sophus::SE3Group<T> TransformT;
    typedef Eigen::Map<TransformT> TransformMapT;
    typedef Eigen::Map<const TransformT> ConstTransformMapT;
    typedef Eigen::Map<RotationT> RotationMapT;
    typedef Eigen::Map<const RotationT> ConstRotationMapT;
    
    typedef Sophus::SO2Group<T> Rotation2DT;
    typedef Sophus::SE2Group<T> Transform2DT;
    typedef Eigen::Map<Transform2DT> Transform2DMapT;
    typedef Eigen::Map<const Transform2DT> ConstTransform2DMapT;
    typedef Eigen::Map<Rotation2DT> Rotation2DMapT;
    typedef Eigen::Map<const Rotation2DT> ConstRotation2DMapT;
    
    typedef typename Sophus::SE3Group<T>::Tangent TangentTransformT;
    typedef Eigen::Map<typename Sophus::SE3Group<T>::Tangent> TangentTransformMapT;
    typedef Eigen::Map<const typename Sophus::SE3Group<T>::Tangent> ConstTangentTransformMapT;
    
    typedef typename Sophus::SE2Group<T>::Tangent TangentTransform2DT;
    typedef Eigen::Map<typename Sophus::SE2Group<T>::Tangent> TangentTransform2DMapT;
    typedef Eigen::Map<const typename Sophus::SE2Group<T>::Tangent> ConstTangentTransform2DMapT;
    
    typedef typename Sophus::SO3Group<T>::Tangent TangentRotationT;
    typedef Eigen::Map<typename Sophus::SO3Group<T>::Tangent> TangentRotationMapT;
    typedef Eigen::Map<const typename Sophus::SO3Group<T>::Tangent> ConstTangentRotationMapT;
    
    typedef typename Sophus::SO2Group<T>::Tangent TangentRotation2DT;
    typedef Eigen::Map<typename Sophus::SO2Group<T>::Tangent> TangentRotation2DMapT;
    typedef Eigen::Map<const typename Sophus::SO2Group<T>::Tangent> ConstTangentRotation2DMapT;
    
    typedef Eigen::Matrix<T,2,1> PixelT;
    typedef Eigen::Map<PixelT> PixelMapT;
    typedef Eigen::Map<const PixelT> ConstPixelMapT;
    
    typedef typename TransformT::Point PointT;
    typedef Eigen::Map<PointT> PointMapT;
    typedef Eigen::Map<const PointT> ConstPointMapT;
    
    typedef typename Transform2DT::Point Point2DT;
    typedef Eigen::Map<Point2DT> Point2DMapT;
    typedef Eigen::Map<const Point2DT> ConstPoint2DMapT;
    
    typedef Eigen::Quaternion<T> QuaternionT;
    typedef Eigen::Map<QuaternionT> QuaternionMapT;
    typedef Eigen::Map<const QuaternionT> ConstQuaternionMapT;
};

template<typename T>
EIGEN_DEVICE_FUNC static inline T getFieldOfView(T focal, T width)
{
    return T(2.0) * atan(width / (T(2.0) * focal));
}

template<typename T>
EIGEN_DEVICE_FUNC static inline T getFocalLength(T fov, T width)
{
    return width / (T(2.0) * tan(fov / T(2.0)));
}

/**
 * Runtime polymorphic interface if one prefers that.
 */
template<typename T>
class CameraInterface
{
public:
    virtual ~CameraInterface() { }
    
    virtual CameraModelType getModelType() const = 0;
    virtual const char* getModelName() const = 0;
    virtual bool pixelValid(T x, T y) const = 0;
    virtual bool pixelValidSquare(T x, T y) const = 0;
    virtual bool pixelValidSquare(const typename ComplexTypes<T>::PixelT& pt) const = 0;
    virtual bool pixelValidCircular(T x, T y) const = 0;
    virtual bool pixelValidCircular(const typename ComplexTypes<T>::PixelT& pt) const = 0;
    virtual typename ComplexTypes<T>::PixelT forward(const typename ComplexTypes<T>::PointT& tmp_pt) const = 0;
    virtual typename ComplexTypes<T>::PointT inverse(T x, T y) const = 0;
    virtual typename ComplexTypes<T>::PixelT forward(const typename ComplexTypes<T>::TransformT& pose, 
                                                     const typename ComplexTypes<T>::PointT& pt) const = 0;
    virtual typename ComplexTypes<T>::PixelT forward(const typename ComplexTypes<T>::RotationT& pose, 
                                                    const typename ComplexTypes<T>::PointT& pt) const = 0;
    virtual typename ComplexTypes<T>::PointT inverse(const typename ComplexTypes<T>::PixelT& pix) const = 0;
    virtual typename ComplexTypes<T>::PointT inverseAtDistance(const typename ComplexTypes<T>::PixelT& pix, T dist) const = 0;
    virtual typename ComplexTypes<T>::PointT inverseAtDistance(T x, T y, T dist) const = 0;
    virtual typename ComplexTypes<T>::PointT inverseAtDistance(const typename ComplexTypes<T>::TransformT& pose, 
                                                               const typename ComplexTypes<T>::PixelT& pix, T dist) const = 0;
    virtual typename ComplexTypes<T>::PointT inverseAtDistance(const typename ComplexTypes<T>::TransformT& pose, 
                                                               T x, T y, T dist) const = 0;
    virtual typename ComplexTypes<T>::PixelT twoFrameProject(const typename ComplexTypes<T>::TransformT& pose1, 
                                                             const typename ComplexTypes<T>::PixelT& pix, T dist, 
                                                             const typename ComplexTypes<T>::TransformT& pose2) const = 0;
    virtual typename ComplexTypes<T>::PixelT twoFrameProject(const typename ComplexTypes<T>::TransformT& pose1, 
                                                             T x, T y, T dist, const typename ComplexTypes<T>::TransformT& pose2) const = 0;
};

/**
 * Functions / methods shared by all the camera models.
 */
template<typename Derived>
class CameraFunctions 
{
public:
    typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;

    // ------------------- statics ---------------------    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT worldToCamera(const typename ComplexTypes<T>::TransformT& pose, const typename ComplexTypes<T>::PointT& pt) 
    {
        return pose.inverse() * pt; // why!
    }

    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT cameraToWorld(const typename ComplexTypes<T>::TransformT& pose, const typename ComplexTypes<T>::PointT& pt) 
    {
        return pose * pt; // why!
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT worldToCamera(const typename ComplexTypes<T>::RotationT& pose, const typename ComplexTypes<T>::PointT& pt) 
    {
        return pose.inverse() * pt; // why!
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT cameraToWorld(const typename ComplexTypes<T>::RotationT& pose, const typename ComplexTypes<T>::PointT& pt) 
    {
        return pose * pt; // why!
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValid(const Derived& ccd, const typename ComplexTypes<T>::PixelT& pt)
    {
        return Derived::template pixelValidSquare<T>(ccd, pt(0), pt(1));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidCircular(const Derived& ccd, const typename ComplexTypes<T>::PixelT& pt)
    {
        return Derived::template pixelValidCircular<T>(ccd, pt(0), pt(1));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidSquare(const Derived& ccd, const typename ComplexTypes<T>::PixelT& pt)
    {
        return Derived::template pixelValidSquare<T>(ccd, pt(0), pt(1));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void resizeViewport(Derived& ccd, const typename ComplexTypes<T>::PixelT& pt)
    {
        Derived::template resizeViewport<T>(ccd, pt(0), pt(1));
    }

    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const Derived& ccd, 
                                                                                          const typename ComplexTypes<T>::TransformT& pose, 
                                                                                          const typename ComplexTypes<T>::PointT& pt)
    {
        return Derived::template forward<T>(ccd, worldToCamera<T>(pose, pt));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const Derived& ccd, 
                                                                                          const typename ComplexTypes<T>::RotationT& pose, 
                                                                                          const typename ComplexTypes<T>::PointT& pt)
    {
        return Derived::template forward<T>(ccd, worldToCamera<T>(pose, pt));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverse(const Derived& ccd, 
                                                                                          const typename ComplexTypes<T>::PixelT& pix)
    {
        return Derived::template inverse<T>(ccd, pix(0), pix(1));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(const Derived& ccd, const typename ComplexTypes<T>::PixelT& pix, T dist)
    {
        return Derived::template inverse<T>(ccd, pix(0), pix(1)) * dist;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(const Derived& ccd, T x, T y, T dist)
    {
        return Derived::template inverse<T>(ccd, x, y) * dist;
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(const Derived& ccd, 
                                                                                                    const typename ComplexTypes<T>::TransformT& pose, 
                                                                                                    const typename ComplexTypes<T>::PixelT& pix, T dist)
    {
        return cameraToWorld<T>(pose, ((Derived::template inverse<T>(ccd,pix)) * dist));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(const Derived& ccd, 
                                                                                                    const typename ComplexTypes<T>::TransformT& pose, 
                                                                                                    T x, T y, T dist)
    {
        return cameraToWorld<T>(pose, ((Derived::template inverse<T>(ccd,x,y)) * dist));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT twoFrameProject(const Derived& ccd, 
                                                                                                  const typename ComplexTypes<T>::TransformT& pose1, 
                                                                                                  const typename ComplexTypes<T>::PixelT& pix, 
                                                                                                  T dist, 
                                                                                                  const typename ComplexTypes<T>::TransformT& pose2)
    {
        return forward<T>(ccd, pose2, inverseAtDistance<T>(ccd, pose1, pix, dist));
    }
    
    template<typename T = Scalar>
    static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT twoFrameProject(const Derived& ccd, 
                                                                                                  const typename ComplexTypes<T>::TransformT& pose1, 
                                                                                                  T x, T y, T dist, 
                                                                                                  const typename ComplexTypes<T>::TransformT& pose2)
    {
        return forward<T>(ccd, pose2, inverseAtDistance<T>(ccd, pose1, x, y, dist));
    }
    
    // ------------------- non statics ---------------------
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void resizeViewport(T x, T y)
    {
        Derived::template resizeViewport<T>(*static_cast<Derived*>(this), x, y);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void resizeViewport(const typename ComplexTypes<T>::PixelT& pt)
    {
        Derived::template resizeViewport<T>(*static_cast<Derived*>(this), pt(0), pt(1));
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValid(T x, T y) const
    {
        return Derived::template pixelValidSquare<T>(*static_cast<const Derived*>(this), x, y);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidSquare(T x, T y) const
    {
        return Derived::template pixelValidSquare<T>(*static_cast<const Derived*>(this), x, y);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidSquare(const typename ComplexTypes<T>::PixelT& pt) const
    {
        return Derived::template pixelValidSquare<T>(*static_cast<const Derived*>(this), pt(0), pt(1));
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidCircular(T x, T y) const
    {
        return Derived::template pixelValidCircular<T>(*static_cast<const Derived*>(this), x, y);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool pixelValidCircular(const typename ComplexTypes<T>::PixelT& pt) const
    {
        return Derived::template pixelValidCircular<T>(*static_cast<const Derived*>(this), pt(0), pt(1));
    }
        
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const typename ComplexTypes<T>::PointT& tmp_pt) const
    {
        return Derived::template forward<T>(*static_cast<const Derived*>(this), tmp_pt);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverse(T x, T y) const
    {
        return Derived::template inverse<T>(*static_cast<const Derived*>(this), x, y);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const typename ComplexTypes<T>::TransformT& pose, const typename ComplexTypes<T>::PointT& pt) const
    {
        return CameraFunctions::forward<T>(*static_cast<const Derived*>(this), pose, pt);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT forward(const typename ComplexTypes<T>::RotationT& pose, const typename ComplexTypes<T>::PointT& pt) const
    {
        return CameraFunctions::forward<T>(*static_cast<const Derived*>(this), pose, pt);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverse(const typename ComplexTypes<T>::PixelT& pix) const
    {
        return CameraFunctions::inverse<T>(*static_cast<const Derived*>(this), pix);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(const typename ComplexTypes<T>::PixelT& pix, T dist) const
    {
        return CameraFunctions::inverseAtDistance<T>(*static_cast<const Derived*>(this), pix, dist); 
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(T x, T y, T dist) const
    {
        return CameraFunctions::inverseAtDistance<T>(*static_cast<const Derived*>(this), x, y, dist);
    }
    
    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(const typename ComplexTypes<T>::TransformT& pose, 
                                                                                const typename ComplexTypes<T>::PixelT& pix, T dist) const
    {
        return CameraFunctions::inverseAtDistance<T>(*static_cast<const Derived*>(this), pose, pix, dist);
    }

    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PointT inverseAtDistance(const typename ComplexTypes<T>::TransformT& pose, T x, T y, T dist) const
    {
        return CameraFunctions::inverseAtDistance<T>(*static_cast<const Derived*>(this), pose, x, y, dist);
    }

    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT twoFrameProject(const typename ComplexTypes<T>::TransformT& pose1, const typename ComplexTypes<T>::PixelT& pix, T dist, 
                                                                                           const typename ComplexTypes<T>::TransformT& pose2) const
    {
        return CameraFunctions::twoFrameProject<T>(*static_cast<const Derived*>(this), pose1, pix, dist, pose2);
    }

    template<typename T = Scalar>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE typename ComplexTypes<T>::PixelT twoFrameProject(const typename ComplexTypes<T>::TransformT& pose1, T x, T y, T dist, 
                                                                                           const typename ComplexTypes<T>::TransformT& pose2) const
    {
        return CameraFunctions::twoFrameProject<T>(*static_cast<const Derived*>(this), pose1, x, y, dist, pose2);
    }
};

/**
 * Wraps templated typed camera model into a polymorphic class.
 */
template <typename ModelT>
class CameraFromCRTP : public ModelT, public CameraInterface<typename ModelT::Scalar>
{
public:
    typedef ModelT Derived;
    typedef typename Derived::Scalar Scalar;
    
    CameraFromCRTP() = default;
    CameraFromCRTP(const Derived& d) : Derived(d) { }
    CameraFromCRTP(const CameraFromCRTP& other) = default;
    
    virtual ~CameraFromCRTP() { }
    
    virtual CameraModelType getModelType() const 
    { 
        return Derived::ModelType; 
    } 
    
    virtual const char* getModelName() const
    {
        return CameraModelToTypeAndName<Derived::ModelType>::Name;
    }
    
    virtual bool pixelValid(Scalar x, Scalar y) const 
    { 
        return Derived::template pixelValid<Scalar>(x,y); 
    }
    
    virtual bool pixelValidSquare(Scalar x, Scalar y) const 
    { 
        return Derived::template pixelValidSquare<Scalar>(x,y); 
    }
    
    virtual bool pixelValidSquare(const typename ComplexTypes<Scalar>::PixelT& pt) const 
    { 
        return Derived::template pixelValidSquare<Scalar>(pt); 
    }
    
    virtual bool pixelValidCircular(Scalar x, Scalar y) const 
    { 
        return Derived::template pixelValidCircular<Scalar>(x,y); 
    }
    
    virtual bool pixelValidCircular(const typename ComplexTypes<Scalar>::PixelT& pt) const 
    { 
        return Derived::template pixelValidCircular<Scalar>(pt); 
    }
    
    virtual typename ComplexTypes<Scalar>::PixelT forward(const typename ComplexTypes<Scalar>::PointT& tmp_pt) const 
    { 
        return Derived::template forward<Scalar>(tmp_pt); 
    }
    
    virtual typename ComplexTypes<Scalar>::PointT inverse(Scalar x, Scalar y) const 
    { 
        return Derived::template inverse<Scalar>(x,y); 
    }
    
    virtual typename ComplexTypes<Scalar>::PixelT forward(const typename ComplexTypes<Scalar>::TransformT& pose, 
                                                          const typename ComplexTypes<Scalar>::PointT& pt) const 
    { 
        return Derived::template forward<Scalar>(pose, pt); 
    }
    
    virtual typename ComplexTypes<Scalar>::PixelT forward(const typename ComplexTypes<Scalar>::RotationT& pose, 
                                                          const typename ComplexTypes<Scalar>::PointT& pt) const 
    { 
        return Derived::template forward<Scalar>(pose, pt); 
    }
    
    virtual typename ComplexTypes<Scalar>::PointT inverse(const typename ComplexTypes<Scalar>::PixelT& pix) const 
    { 
        return Derived::template inverse<Scalar>(pix); 
    }
    
    virtual typename ComplexTypes<Scalar>::PointT inverseAtDistance(const typename ComplexTypes<Scalar>::PixelT& pix, 
                                                                    Scalar dist) const 
    { 
        return Derived::template inverseAtDistance<Scalar>(pix,dist); 
    }
    
    virtual typename ComplexTypes<Scalar>::PointT inverseAtDistance(Scalar x, Scalar y, Scalar dist) const 
    { 
        return Derived::template inverseAtDistance<Scalar>(x,y,dist); 
    }
    
    virtual typename ComplexTypes<Scalar>::PointT inverseAtDistance(const typename ComplexTypes<Scalar>::TransformT& pose, 
                                                                    const typename ComplexTypes<Scalar>::PixelT& pix, 
                                                                    Scalar dist) const 
    { 
        return Derived::template inverseAtDistance<Scalar>(pose,pix,dist); 
    }
    
    virtual typename ComplexTypes<Scalar>::PointT inverseAtDistance(const typename ComplexTypes<Scalar>::TransformT& pose, 
                                                                    Scalar x, Scalar y, Scalar dist) const 
    { 
        return Derived::template inverseAtDistance<Scalar>(pose,x,y,dist); 
    }
    
    virtual typename ComplexTypes<Scalar>::PixelT twoFrameProject(const typename ComplexTypes<Scalar>::TransformT& pose1, 
                                                                  const typename ComplexTypes<Scalar>::PixelT& pix, 
                                                                  Scalar dist, 
                                                                  const typename ComplexTypes<Scalar>::TransformT& pose2) const 
    { 
        return Derived::template twoFrameProject<Scalar>(pose1,pix,dist,pose2); 
    }
    
    virtual typename ComplexTypes<Scalar>::PixelT twoFrameProject(const typename ComplexTypes<Scalar>::TransformT& pose1, 
                                                                  Scalar x, Scalar y, Scalar dist, 
                                                                  const typename ComplexTypes<Scalar>::TransformT& pose2) const 
    { 
        return Derived::template twoFrameProject<Scalar>(pose1,x,y,dist,pose2); 
    }
};

}

#endif // CAMERA_MODEL_HELPERS_HPP
