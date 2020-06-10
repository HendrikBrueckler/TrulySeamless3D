#include <stdlib.h>
#include <time.h>

#include "TrulySeamless3D/trulyseamless.h"

int main(int argc, char** argv)
{
    if (argc > 1)
    {
        std::string fileName(argv[1]);
        TS3D::TrulySeamless3D sanitizer(fileName);
        bool keepOriginalTransitions = false;
        if (argc == 3 && std::string(argv[2]) == "--preserve-cut-graph")
            keepOriginalTransitions = true;
        if (!sanitizer.sanitize(0.0, keepOriginalTransitions))
        {
            std::cout << "Could NOT sanitize input, check your input (or file a bug report)" << std::endl;
            return 1;
        }
        sanitizer.writeToFile(fileName + "-seamless.hexex");
        return 0;
    }
    else
    {
        std::cout << "Usage: TrulySeamless3D <Path to tet mesh> [--preserve-cut-graph]" << std::endl;
        return -1;
    }
}
