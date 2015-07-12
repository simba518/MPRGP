#include "mainConfig.h"
#ifdef MAIN_FILE_FORMAT

#include "FEMMeshFormat.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

void main(int argc,char** argv)
{
    if(argc != 3) {
        NOTIFY_MSG("Usage: exe [input] [output]");
        return;
    }

    boost::filesystem::path in(argv[1]);
    boost::filesystem::path out(argv[2]);

    if(in.extension().string() == ".mesh" ||
            in.extension().string() == ".MESH")
        FEMMeshFormat::meshToABQ(in.string(),out.string());

    if(in.extension().string() == ".vtk" ||
            in.extension().string() == ".VTK")
        FEMMeshFormat::VTKToABQ(in.string(),out.string());
}

#endif
