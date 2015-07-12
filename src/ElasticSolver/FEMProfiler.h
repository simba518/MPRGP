#ifndef PROFILER_H
#define PROFILER_H

#include "FEMSystem.h"
#include <boost/filesystem/fstream.hpp>
#include <time.h>

PRJ_BEGIN

class FEMProfiler
{
public:
    typedef Eigen::Matrix<scalarD,-1,1> Vec;
    FEMProfiler(const std::string& path,FEMSolver* sol):_os(path),_sol(sol) {
        _os << "EK,EP,EK+EP,TIME(ms)" << std::endl;
    }
    void beginFrame() {
        _time=clock();
    }
    void endFrame() {
        scalarD TMS=clock()-_time;
        scalarD EK,EP;
        _sol->getSystemEnergy(EK,EP);
        _os<< EK << "," << EP << "," << EK+EP << "," << TMS << std::endl;
    }
protected:
    boost::filesystem::ofstream _os;
    std::clock_t _time;
    FEMSolver* _sol;
};

PRJ_END

#endif