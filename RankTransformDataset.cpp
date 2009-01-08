#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace bll = boost::lambda;

class RankTransformer
{
public:
  RankTransformer(const std::string& aMatrixDir, const std::string& aOutputfile,
                  bool aQuantileNormalisation = false, bool aUseInverse = false)
    : mMatrixDir(aMatrixDir), mOutputFile(NULL), mData(NULL), mBuf(NULL),
      mRanks(NULL), mInvRanks(NULL), mRankAvgs(NULL), mRankCounts(NULL),
      mQuantileNormalisation(aQuantileNormalisation), mUseInverse(aUseInverse)
  {
    mOutputFile = fopen(aOutputfile.c_str(), "w");

    fs::path data(mMatrixDir);
    if (mUseInverse)
      data /= "inverse_data";
    else
      data /= "data";

    mData = fopen(data.string().c_str(), "r");

    fs::path genes(mMatrixDir);
    // Ugly hack: if we are using the inverted data, we simply swap out the
    // list of genes for the list of arrays, so that nGenes is actually the
    // number of arrays. This means that the normalisation occurs as normal,
    // except for each gene across arrays instead of for each array across
    // genes.
    if (mUseInverse)
      genes /= "arrays";
    else
      genes /= "genes";
    nGenes = 0;
    std::ifstream gs(genes.string().c_str());
    while (gs.good())
    {
      std::string l;
      std::getline(gs, l);
      if (!gs.good())
        break;

      nGenes++;
    }

    mBuf = new double[nGenes];
    mRanks = new double[nGenes];
    if (mQuantileNormalisation)
    {
      mRankAvgs = new double[nGenes];
      mRankCounts = new uint32_t[nGenes];
      memset(mRankAvgs, 0, sizeof(double) * nGenes);
      memset(mRankCounts, 0, sizeof(uint32_t) * nGenes);
    }
    mInvRanks = new uint32_t[nGenes];

    processAllData();
  }

  ~RankTransformer()
  {
    if (mOutputFile != NULL)
      fclose(mOutputFile);
    if (mData != NULL)
      fclose(mData);
    if (mBuf != NULL)
      delete [] mBuf;
    if (mInvRanks != NULL)
      delete [] mInvRanks;
    if (mRankCounts != NULL)
      delete [] mRankCounts;
    if (mRankAvgs != NULL)
      delete [] mRankAvgs;
  }

private:
  fs::path mMatrixDir;
  FILE* mOutputFile, * mData;
  double* mBuf, * mRanks;
  uint32_t nGenes;
  uint32_t* mInvRanks;
  double * mRankAvgs;
  uint32_t * mRankCounts;
  bool mQuantileNormalisation;
  bool mUseInverse;

  void
  processAllData()
  {
    if (mQuantileNormalisation)
    {
      while (fread(mBuf, sizeof(double), nGenes, mData) == nGenes)
      {
        processArray(true);
      }
      fseek(mData, 0, SEEK_SET);
      for (uint32_t i = 0; i < nGenes; i++)
        mRankAvgs[i] /= mRankCounts[i];
    }
    while (fread(mBuf, sizeof(double), nGenes, mData) == nGenes)
    {
      processArray();
    }
  }

  void
  processArray(bool aAverageMode = false)
  {
    uint32_t nNans = 0;
    for (uint32_t i = 0; i < nGenes; i++)
    {
      mInvRanks[i] = i;
      nNans += !!finite(mBuf[i]);
    }

    uint32_t nNotNans = nGenes - nNans;
    double rankInflationFactor = (nGenes + 0.0) / nNotNans;

    // Sort indices by value, with NaNs at the top.
    std::sort(mInvRanks, mInvRanks + nGenes,
              (bll::var(mBuf)[bll::_1] < bll::var(mBuf)[bll::_2]) ||
              (bll::bind(finite, bll::var(mBuf)[bll::_1]) &&
               !bll::bind(finite, bll::var(mBuf)[bll::_2])));
    
    uint32_t i;
    if (aAverageMode)
    {
      for (i = 0; i < nGenes && finite(mBuf[mInvRanks[i]]); i++)
      {
        mRankAvgs[i] += mBuf[mInvRanks[i]];
        mRankCounts[i]++;
      }
    }
    else
    {
      if (mQuantileNormalisation)
      {
        for (i = 0; i < nGenes && finite(mBuf[mInvRanks[i]]); i++)
          // We could put code in here to deal with tied ranks by putting in median
          // ranks, but I doubt it would make enough difference to justify it.
          mRanks[mInvRanks[i]] = mRankAvgs[i];
      }
      else
      {
        for (i = 0; i < nGenes && finite(mBuf[mInvRanks[i]]); i++)
          // We could put code in here to deal with tied ranks by putting in median
          // ranks, but I doubt it would make enough difference to justify it.
          mRanks[mInvRanks[i]] = i * rankInflationFactor;
      }

      for (; i < nGenes; i++)
        mRanks[mInvRanks[i]] = std::numeric_limits<double>::quiet_NaN();

      fwrite(mRanks, sizeof(double), nGenes, mOutputFile);
    }
  }
};

int
main(int argc, char** argv)
{
  std::string matrixdir, outputfile;
  po::options_description desc;

  desc.add_options()
    ("matrixdir", po::value<std::string>(&matrixdir), "The directory to read the data from")
    ("use_inverse", "Use the inverted data-set instead of the original")
    ("output", po::value<std::string>(&outputfile), "The file to write the output into")
    ("qnorm", "If specified, causes quantile normalisation to be applied to the data")
    ;

  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("matrixdir"))
      wrong = "matrixdir";
    else if (!vm.count("output"))
      wrong = "output";
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
    std::cerr << "Invalid matrix directory path supplied" << std::endl;
    return 1;
  }

  RankTransformer rt(matrixdir, outputfile, vm.count("qnorm") != 0, vm.count("use_inverse") != 0);
}
