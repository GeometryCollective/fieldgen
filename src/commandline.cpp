#include <string>
#include <iostream>
using namespace std;

#include "DenseMatrix.h"
#include "Mesh.h"
using namespace DDG;

// Command line parsing code courtesy Rohan Sawhney!

void printUsage(const std::string& programName)
{
   std::cout << "usage: "
             << programName << " "
             << "OBJ_INPUT_PATH "
             << "OBJ_OUTPUT_PATH "
             << "[--degree=n] "
             << "[--alignToCurvature] "
             << "[--alignToBoundary] "
             << "[--alignToGivenField] "
             << "[--s=S] "
             << "[--t=T]"
             << "[--alignmentMagnitude=alignmentMagnitude]"
             << std::endl;
}

bool doesArgExist(const std::string& arg, const std::string& searchStr)
{
   return arg.find(searchStr) != std::string::npos;
}

bool parseArg(const std::string& arg, const std::string& searchStr, std::string& value)
{
   if (doesArgExist(arg, searchStr)) {
      value = arg.substr(arg.find_first_of(searchStr[searchStr.size()-1]) + 1);
      return true;
   }

   return false;
}

void parseArgs(int argc, char *argv[], std::string& inputPath, std::string& outputPath,
      int& degree, bool& alignToCurvature, bool& alignToBoundary, bool& alignToGivenField,
      double& s, double& t, double& alignmentMagnitude)
{
   if (argc < 3) {
      // input and/or output path not specified
      printUsage(argv[0]);
      exit(EXIT_FAILURE);

   } else {
      // parse arguments
      inputPath = argv[1];
      outputPath = argv[2];
      std::string degreeStr;
      std::string sStr, tStr, alignmentMagnitudeStr;

      for (int i = 3; i < argc; i++) {
         if (parseArg(argv[i], "--degree=", degreeStr)) degree = atoi(degreeStr.c_str());
         if (doesArgExist(argv[i], "--alignToCurvature")) alignToCurvature = true;
         if (doesArgExist(argv[i], "--alignToBoundary")) alignToBoundary = true;
         if (doesArgExist(argv[i], "--alignToGivenField")) alignToGivenField = true;
         if (parseArg(argv[i], "--s=", sStr)) s = std::atof(sStr.c_str());
         if (parseArg(argv[i], "--t=", tStr)) t = std::atof(tStr.c_str());
         if (parseArg(argv[i], "--alignmentMagnitude=", alignmentMagnitudeStr))
            alignmentMagnitude = std::atof(alignmentMagnitudeStr.c_str());
      }
   }

   // aligning to boundary takes precedence over aligning to curvature
   if (alignToBoundary) {
      alignToCurvature = false;
   }
}

int main( int argc, char** argv )
{
   // parse command line options
   std::string inputPath = "";
   std::string outputPath = "";
   int degree = 1;
   bool alignToCurvature = false;
   bool alignToBoundary = false;
   bool alignToGivenField = false; // the parameter to be passed when to turn on given alignment
   double s = 0.;
   double t = 0.;
   double alignmentMagnitude = 1.0; // not being used right now
   parseArgs( argc, argv, inputPath, outputPath, degree, alignToCurvature, alignToBoundary, alignToGivenField, s, t, alignmentMagnitude );

   Mesh mesh;

   cout << "Reading mesh from " << inputPath << "..." << endl;
   mesh.read( inputPath );

   cout << "Computing field..." << endl;
   mesh.InitKVecDirData();
   mesh.clearSingularities();
   if( alignToCurvature )
   {
      mesh.SmoothestCurvatureAlignment( degree, s, t, true );
   } 
   else if( alignToGivenField )
   {
      mesh.SmoothestGivenVectorAlignment( degree, s, t, alignmentMagnitude, true );
   }
   else if( alignToBoundary )
   {
      mesh.ComputeSmoothestFixedBoundary( degree, s, true );
   }
   else
   {
      mesh.ComputeSmoothest( degree, s, true );
   }

   cout << "Writing solution to " << outputPath << "..." << endl;
   mesh.write( outputPath, degree );

   return 0;
}

