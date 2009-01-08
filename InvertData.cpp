#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class RuntimeException
{
public:
  RuntimeException(const char* aWhy)
    : mWhy(aWhy)
  {
  }

  const std::string what()
  {
    return mWhy;
  }

private:
  const char* mWhy;
};

class DataInverter
{
public:
  DataInverter(const std::string& aMatrixDir)
    : mMatrixDir(aMatrixDir), mnArrays(0), mnGenes(0)
  {
    fs::path arrayList(mMatrixDir);
    arrayList /= "arrays";
    {
      std::ifstream as(arrayList.string().c_str());
      while (as.good())
      {
        std::string l;
        std::getline(as, l);
        if (!as.good())
          break;
        mnArrays++;
      }
    }

    fs::path geneList(mMatrixDir);
    geneList /= "genes";
    {
      std::ifstream gs(geneList.string().c_str());
      while (gs.good())
      {
        std::string l;
        std::getline(gs, l);
        if (!gs.good())
          break;
        mnGenes++;
      }
    }

    fs::path data(mMatrixDir);
    data /= "data";

    fs::path invdata(mMatrixDir);
    invdata /= "inverse_data";
    FILE* dataf = fopen(data.string().c_str(), "r");
    FILE* invdataf = fopen(invdata.string().c_str(), "w");

    double* bigbuf = new double[mnGenes * kConcurrentRows];
    double* smallbuf = new double[kConcurrentRows];
    uint32_t row0 = 0;
    while (row0 < mnArrays)
    {
      uint32_t rownext = row0 + kConcurrentRows;
      if (rownext > mnArrays)
        rownext = mnArrays;
      if (fread(bigbuf, (rownext - row0) * mnGenes, sizeof(double), dataf) != sizeof(double))
      {
        delete [] bigbuf;
        delete [] smallbuf;
        throw RuntimeException("data file is truncated.");
      }

      fseek(invdataf, row0 * sizeof(double), SEEK_SET);
 
      for (uint32_t col = 0; col < mnGenes; col++)
      {
        const double* p = bigbuf + col;
        for (uint32_t i = 0; i < (rownext - row0); i++)
          smallbuf[i] = p[mnGenes * i];
        fwrite(smallbuf, (rownext - row0), sizeof(double), invdataf);
        fseek(invdataf, (mnArrays - rownext + row0) * sizeof(double), SEEK_CUR);
      }

      row0 = rownext;
    }

    delete [] bigbuf;
    delete [] smallbuf;
  }

private:
  static const uint32_t kConcurrentRows = 3000;
  std::string mMatrixDir;
  uint32_t mnArrays, mnGenes;
};

int
main(int argc, char**argv)
{
  std::string matrixdir;
  po::options_description desc;

  desc.add_options()
    ("matrixdir", po::value<std::string>(&matrixdir),
     "write matrix into directory")
    ("help", "produce help message")
    ;

  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("matrixdir"))
      wrong = "matrixdir";
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;
  if (vm.count("help") || wrong != "")
  {
    std::cout << desc << std::endl;
    return 1;
  }
  
  if (!fs::is_directory(matrixdir))
  {
    std::cerr << "Matrix 'directory' is not a directory." << std::endl;
    return 1;
  }

  try
  {
    DataInverter di(matrixdir);
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}
