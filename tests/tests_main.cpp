/**
 * 
 * Generalized Camera Models.
 * GTest Main Launcher.
 * 
 * Copyright (c) Robert Lukierski 2015. All rights reserved.
 * Author: Robert Lukierski.
 * 
 */

// testing framework & libraries
#include <gtest/gtest.h>

// google logger
#include <glog/logging.h>
#include <gflags/gflags.h>

extern "C" int main(int argc, char** argv)
{
    google::InitGoogleLogging(argv[0]);
    google::LogToStderr(); // no files in /tmp
    google::SetStderrLogging(google::INFO);
    
    ::testing::InitGoogleTest(&argc, argv);
    int rc = RUN_ALL_TESTS();
    
    google::ShutdownGoogleLogging(); 
    gflags::ShutDownCommandLineFlags();
    
    return rc;
}
