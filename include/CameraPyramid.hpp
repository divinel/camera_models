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
 * Pyramid for camera models to rescale coefficients accordingly.
 * ****************************************************************************
 */

#ifndef CAMERA_PYRAMID_HPP
#define CAMERA_PYRAMID_HPP

#include <vector>
#include <CameraModelHelpers.hpp>

namespace camera
{

// Power of two pyramid.
template<typename MODEL_T, std::size_t Levels>
class CameraPyramid
{
public:
    static const std::size_t LevelCount = Levels;
    typedef MODEL_T CameraModelT;
    typedef typename CameraModelT::Scalar Scalar;
    
    EIGEN_DEVICE_FUNC inline CameraPyramid()
    {
        
    }
    
    inline CameraPyramid(const CameraModelT& initial)
    {
        models[0] = initial;
        pyramidDown();
    }
    
    EIGEN_DEVICE_FUNC inline ~CameraPyramid()
    {
        
    }
    
    EIGEN_DEVICE_FUNC inline CameraPyramid(const CameraPyramid<CameraModelT,LevelCount>& pyramid)
    {
        for(std::size_t l = 0 ; l < LevelCount ; ++l) 
        {
            models[l] = pyramid.models[l];
        }
    }
    
    EIGEN_DEVICE_FUNC inline CameraPyramid(CameraPyramid<CameraModelT,LevelCount>&& pyramid) 
    {
        for(std::size_t l = 0 ; l < LevelCount ; ++l) 
        {
            models[l] = std::move(pyramid.models[l]);
        }
    }
    
    EIGEN_DEVICE_FUNC inline CameraPyramid<CameraModelT,LevelCount>& operator=(const CameraPyramid<CameraModelT,LevelCount>& pyramid)
    {
        for(std::size_t l = 0 ; l < LevelCount ; ++l) 
        {
            models[l] = pyramid.models[l];
        }

        return *this;
    }
    
    EIGEN_DEVICE_FUNC inline CameraPyramid<CameraModelT,LevelCount>& operator=(CameraPyramid<CameraModelT,LevelCount>&& pyramid)
    {
        for(std::size_t l = 0 ; l < LevelCount ; ++l) 
        {
            models[l] = std::move(pyramid.models[l]);
        }
        
        return *this;
    }
    
    inline void pyramidDown()
    {
        for(std::size_t l = 1; l < LevelCount ; ++l ) 
        {
            CameraModelT prev = models[l-1];
            prev.resizeViewport(prev.width()/Scalar(2.0),prev.height()/Scalar(2.0));
            models[l] = prev;
        }
    }

    EIGEN_DEVICE_FUNC inline CameraModelT& operator[](std::size_t i)
    {
        assert(i < LevelCount);
        return models[i];
    }

    EIGEN_DEVICE_FUNC inline const CameraModelT& operator[](std::size_t i) const
    {
        assert(i < LevelCount);
        return models[i];
    }
    
    EIGEN_DEVICE_FUNC inline CameraModelT& operator()(std::size_t i)
    {
        assert(i < LevelCount);
        return models[i];
    }
    
    EIGEN_DEVICE_FUNC inline const CameraModelT& operator()(std::size_t i) const
    {
        assert(i < LevelCount);
        return models[i];
    }

    template<std::size_t SubLevels>
    EIGEN_DEVICE_FUNC inline CameraPyramid<CameraModelT,SubLevels> subPyramid(std::size_t startLevel)
    {
        assert(startLevel + SubLevels < LevelCount);
        
        CameraPyramid<CameraModelT,SubLevels> pyr;

        for(std::size_t l = 0 ; l < SubLevels; ++l) 
        {
            pyr.models[l] = models[startLevel+l];
        }

        return pyr;
    }
protected:
    CameraModelT models[LevelCount];
};

// Power of two pyramid, runtime number of levels.
template<typename MODEL_T>
class RuntimeCameraPyramid
{
public:
    typedef MODEL_T CameraModelT;
    typedef typename CameraModelT::Scalar Scalar;
    
    inline RuntimeCameraPyramid() = delete;
    
    inline RuntimeCameraPyramid(std::size_t Levels) : models(Levels)
    {
        
    }
    
    inline RuntimeCameraPyramid(std::size_t Levels, const CameraModelT& initial) : models(Levels)
    {
        models[0] = initial;
        pyramidDown();
    }
    
    inline ~RuntimeCameraPyramid()
    {
        
    }
    
    inline RuntimeCameraPyramid(const RuntimeCameraPyramid<CameraModelT>& pyramid)
    {
        models.resize(pyramid.getLevelCount());
        
        for(std::size_t l = 0 ; l < getLevelCount() ; ++l) 
        {
            models[l] = pyramid.models[l];
        }
    }
    
    inline RuntimeCameraPyramid(RuntimeCameraPyramid<CameraModelT>&& pyramid) 
    {
        models.resize(pyramid.getLevelCount());
        
        for(std::size_t l = 0 ; l < getLevelCount() ; ++l) 
        {
            models[l] = std::move(pyramid.models[l]);
        }
    }
    
    inline RuntimeCameraPyramid<CameraModelT>& operator=(const RuntimeCameraPyramid<CameraModelT>& pyramid)
    {
        models.resize(pyramid.getLevelCount());
        
        for(std::size_t l = 0 ; l < getLevelCount() ; ++l) 
        {
            models[l] = pyramid.models[l];
        }

        return *this;
    }
    
    inline RuntimeCameraPyramid<CameraModelT>& operator=(RuntimeCameraPyramid<CameraModelT>&& pyramid)
    {
        models.resize(pyramid.getLevelCount());
        
        for(std::size_t l = 0 ; l < getLevelCount() ; ++l) 
        {
            models[l] = std::move(pyramid.models[l]);
        }
        
        return *this;
    }
    
    inline void swap(RuntimeCameraPyramid<CameraModelT>& pyramid)
    {
        for(std::size_t l = 0 ; l < getLevelCount() ; ++l) 
        {
            models[l].swap(pyramid.models[l]);
        }
    }
    
    inline void pyramidDown()
    {
        for(std::size_t l = 1; l < getLevelCount() ; ++l ) 
        {
            CameraModelT prev = models[l-1];
            prev.resizeViewport(prev.width()/Scalar(2.0),prev.height()/Scalar(2.0));
            models[l] = prev;
        }
    }

    inline CameraModelT& operator[](std::size_t i)
    {
        assert(i < getLevelCount());
        return models[i];
    }

    inline const CameraModelT& operator[](std::size_t i) const
    {
        assert(i < getLevelCount());
        return models[i];
    }
    
    inline CameraModelT& operator()(std::size_t i)
    {
        assert(i < getLevelCount());
        return models[i];
    }
    
    inline const CameraModelT& operator()(std::size_t i) const
    {
        assert(i < getLevelCount());
        return models[i];
    }

    inline RuntimeCameraPyramid<CameraModelT> subPyramid(std::size_t startLevel, std::size_t SubLevels)
    {
        assert(startLevel + SubLevels < getLevelCount());
        
        RuntimeCameraPyramid<CameraModelT> pyr;

        for(std::size_t l = 0 ; l < SubLevels; ++l) 
        {
            pyr.models[l] = models[startLevel+l];
        }

        return pyr;
    }
    
    inline std::size_t getLevelCount() const { return models.size(); }
protected:
    std::vector<CameraModelT> models;
};

}

#endif // CAMERA_PYRAMID_HPP
