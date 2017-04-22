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
 * Projection matrices to be used with OpenGL.
 * ****************************************************************************
 */

#ifndef CAMERA_PROJECTION_MATRICES_HPP
#define CAMERA_PROJECTION_MATRICES_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// http://www.songho.ca/opengl/gl_projectionmatrix.html

namespace camera
{
    
template<typename T>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<T,4,4> getOrtographicProjection(T l, T r, T b, T t, T n, T f)
{
    Eigen::Matrix<T,4,4> m;
    
    m << T(2.0)/(r-l) , T(0.0)       , T(0.0)        , -(r+l)/(r-l),
         T(0.0)       , T(2.0)/(t-b) , T(0.0)        , -(t+b)/(t-b),
         T(0.0)       , T(0.0)       , T(-2.0)/(f-n) , -(f+n)/(f-n),
         T(0.0)       , T(0.0)       , T(0.0)        , T(1.0);

    return m;
}

// Camera Axis:
//   X - Right, Y - Up, Z - Back
// Image Origin:
//   Bottom Left
// Caution: Principal point defined with respect to image origin (0,0) at
//          top left of top-left pixel (not center, and in different frame
//          of reference to projection function image)
template<typename CameraModel>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<typename CameraModel::Scalar,4,4> getPerspectiveProjectionRUBBottomLeft(const CameraModel& cam, typename CameraModel::Scalar N, typename CameraModel::Scalar F)
{
    typedef typename CameraModel::Scalar Scalar;
    
    const Scalar L = +(cam.u0()) * N / -cam.fx();
    const Scalar T = +(cam.v0()) * N / cam.fy();
    const Scalar R = -(cam.width() - cam.u0()) * N / -cam.fx();
    const Scalar B = -(cam.height() - cam.v0()) * N / cam.fy();
 
    Eigen::Matrix<typename CameraModel::Scalar,4,4> ret;
    
    ret << 
    Scalar(2.0)*N/(R-L) , Scalar(0.0)         , (R+L)/(R-L) , Scalar(0.0),
    Scalar(0.0)         , Scalar(2.0)*N/(T-B) , (T+B)/(T-B) , Scalar(0.0),
    Scalar(0.0)         , Scalar(0.0)         ,-(F+N)/(F-N) ,-(Scalar(2.0)*F*N)/(F-N),
    Scalar(0.0)         , Scalar(0.0)         ,-Scalar(1.0) , Scalar(0.0);
    
    return ret;
}

// Camera Axis:
//   X - Right, Y - Down, Z - Forward
// Image Origin:
//   Top Left
// Pricipal point specified with image origin (0,0) at top left of top-left pixel (not center)
template<typename CameraModel>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<typename CameraModel::Scalar,4,4> getPerspectiveProjectionRDFTopLeft(const CameraModel& cam, typename CameraModel::Scalar N, typename CameraModel::Scalar F)
{
    typedef typename CameraModel::Scalar Scalar;
    
    const Scalar L = -(cam.u0()) * N / cam.fx();
    const Scalar T = -(cam.v0()) * N / cam.fy();
    const Scalar R = +(cam.width() - cam.u0()) * N / cam.fx();
    const Scalar B = +(cam.height() - cam.v0()) * N / cam.fy();
        
    Eigen::Matrix<typename CameraModel::Scalar,4,4> ret;
    
    ret << 
    Scalar(2.0)*N/(R-L) , Scalar(0.0)         , (R+L)/(L-R) , Scalar(0.0),
    Scalar(0.0)         , Scalar(2.0)*N/(T-B) , (T+B)/(B-T) , Scalar(0.0),
    Scalar(0.0)         , Scalar(0.0)         , (F+N)/(F-N) ,(Scalar(2.0)*F*N)/(N-F),
    Scalar(0.0)         , Scalar(0.0)         , Scalar(1.0) , Scalar(0.0);
    
    return ret;
}

// Camera Axis:
//   X - Right, Y - Down, Z - Forward
// Image Origin:
//   Bottom Left
// Pricipal point specified with image origin (0,0) at top left of top-left pixel (not center)
template<typename CameraModel>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<typename CameraModel::Scalar,4,4> getPerspectiveProjectionRDFBottomLeft(const CameraModel& cam, typename CameraModel::Scalar N, typename CameraModel::Scalar F)
{
    typedef typename CameraModel::Scalar Scalar;
    
    const Scalar L = -(cam.u0()) * N / cam.fx();
    const Scalar T = +(cam.height() - cam.v0()) * N / cam.fy();
    const Scalar R = +(cam.width() - cam.u0()) * N / cam.fx();
    const Scalar B = -(cam.v0()) * N / cam.fy();
    
    Eigen::Matrix<typename CameraModel::Scalar,4,4> ret;
    
    ret << 
    Scalar(2.0)*N/(R-L) , Scalar(0.0)         , (R+L)/(L-R) , Scalar(0.0),
    Scalar(0.0)         , Scalar(2.0)*N/(T-B) , (T+B)/(B-T) , Scalar(0.0),
    Scalar(0.0)         , Scalar(0.0)         , (F+N)/(F-N) ,(Scalar(2.0)*F*N)/(N-F),
    Scalar(0.0)         , Scalar(0.0)         , Scalar(1.0) , Scalar(0.0);
    
    return ret;
}

template<typename CameraModel>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<typename CameraModel::Scalar,4,4> getPerspectiveProjection(const CameraModel& cam, typename CameraModel::Scalar N, typename CameraModel::Scalar F)
{
    return getPerspectiveProjectionRDFTopLeft(cam,N,F);
}

// NOTE element (2,3) may be wrong

template<typename T>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<T,4,4> lookAtRUB(const Eigen::Matrix<T,3,1>& eye, const Eigen::Matrix<T,3,1>& center, const Eigen::Matrix<T,3,1>& up)
{
    Eigen::Matrix<T,4,4> view_matrix = Eigen::Matrix<T,4,4>::Zero();
    Eigen::Matrix<T,3,3> R;
    R.col(2) = (eye-center).normalized();
    R.col(0) = up.cross(R.col(2)).normalized();
    R.col(1) = R.col(2).cross(R.col(0));
    view_matrix.template topLeftCorner<3,3>() = R.transpose();
    view_matrix.template topRightCorner<3,1>() = -R.transpose() * center;
    view_matrix.row(3) << T(0.0), T(0.0), T(0.0), T(1.0);
    return view_matrix;
}

template<typename T>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<T,4,4> lookAtRDF(const Eigen::Matrix<T,3,1>& eye, const Eigen::Matrix<T,3,1>& center, const Eigen::Matrix<T,3,1>& up)
{
    Eigen::Matrix<T,4,4> view_matrix = Eigen::Matrix<T,4,4>::Zero();
    Eigen::Matrix<T,3,3> R;
    R.col(2) = (center-eye).normalized();
    R.col(0) = R.col(2).cross(up).normalized();
    R.col(1) = R.col(2).cross(R.col(0));
    view_matrix.template topLeftCorner<3,3>() = R.transpose();
    view_matrix.template topRightCorner<3,1>() = -R.transpose() * center;
    view_matrix.row(3) << T(0.0), T(0.0), T(0.0), T(1.0);
    return view_matrix;
}

template<typename T>
EIGEN_DEVICE_FUNC inline Eigen::Matrix<T,4,4> lookAt(const Eigen::Matrix<T,3,1>& eye, const Eigen::Matrix<T,3,1>& center, const Eigen::Matrix<T,3,1>& up)
{
    return lookAtRDF(eye, center, up);
}

}

#endif // CAMERA_PROJECTION_MATRICES_HPP