#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <list>
#include <boost/tokenizer.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class SOFT2Matrix
{
public:
  SOFT2Matrix(const std::string& aSOFTFile, const std::string& aOutdir)
    : mOutdir(aOutdir), mSOFTFile(aSOFTFile.c_str()), mnSamples(0),
      mNextIndex(0)
  {
    fs::path arrayList(mOutdir);
    arrayList /= "arrays";

    mArrayList = new std::ofstream(arrayList.string().c_str());

    fs::path geneList(mOutdir);
    geneList /= "genes";

    mGeneList = new std::ofstream(geneList.string().c_str());
  }

  ~SOFT2Matrix()
  {
    delete mArrayList;
    delete mGeneList;
  }

  void
  process()
  {
    std::string l;
    processLine = &SOFT2Matrix::processPlatformIntro;

    while (!mSOFTFile.bad())
    {
      std::getline(mSOFTFile, l);
      (this->*processLine)(l);
    }
  }

private:
  fs::path mOutdir;
  std::ifstream mSOFTFile;
  std::ofstream *mArrayList, *mGeneList;
  void (SOFT2Matrix::* processLine)(const std::string& aLine);
  uint32_t mnSamples;
  std::list<std::string> mSampleIds;
  std::list<std::string>::iterator mNextId;

  void
  processPlatformIntro(const std::string& aLine)
  {
    if (aLine == "!platform_table_begin")
    {
      processLine = &SOFT2Matrix::processPlatformHeader;
      return;
    }

    if (aLine.substr(0, 22) == "!Platform_sample_id = ")
    {
      std::string sampleId(aLine.substr(22));
      mSampleIds.push_back(sampleId);
      (*mArrayList) << sampleId << std::endl;
      mnSamples++;
      if (mnSamples == 1)
        mNextId = mSampleIds.begin();
    }
  }

  void
  processPlatformHeader(const std::string& aLine)
  {
    processLine = &SOFT2Matrix::processPlatformTable;

    boost::char_separator<char> tdv("\t", "", boost::keep_empty_tokens);
    typedef boost::tokenizer<boost::char_separator<char> > tok_t;
    tok_t tok(aLine, tdv);

    uint32_t n = 0;

    for (tok_t::iterator i = tok.begin(); i != tok.end(); i++, n++)
    {
      if ((*i) == "ID")
        mIdIndex = n;
      else if ((*i) == "Gene Symbol")
        mGeneSymbolIndex = n;
    }
  }

  uint32_t mIdIndex, mGeneSymbolIndex, mValueIndex;
  uint32_t mNextIndex;

  void
  processPlatformTable(const std::string& aLine)
  {
    if (aLine == "!platform_table_end")
    {
      processLine = &SOFT2Matrix::processSampleIntro;
      return;
    }

    boost::char_separator<char> tdv("\t", "", boost::keep_empty_tokens);
    typedef boost::tokenizer<boost::char_separator<char> > tok_t;
    tok_t tok(aLine, tdv);

    uint32_t n = 0;

    std::string id, symbol;
    for (tok_t::iterator i = tok.begin(); i != tok.end(); i++, n++)
    {
      if (n == mIdIndex)
        id = *i;
      else if (n == mGeneSymbolIndex)
        symbol = *i;
    }

    mGeneIndexById.insert(std::pair<std::string, uint32_t>(id, mNextIndex));
    (*mGeneList) << symbol << std::endl;
    mNextIndex++;
  }

  std::map<std::string, uint32_t> mGeneIndexById;

  void
  processSampleIntro(const std::string& aLine)
  {
    if (aLine.substr(0, 10) == "^SAMPLE = ")
    {
      std::string sampId = aLine.substr(10);
      if (sampId != *mNextId)
        std::cout << "Sample ID mismatch: expected "
                  << *mNextId << " got " << sampId
                  << std::endl;
      mNextId++;
      return;
    }
    if (aLine == "!sample_table_begin")
      processLine = &SOFT2Matrix::processSampleHeader;
  }

  void
  processSampleHeader(const std::string& aLine)
  {
    processLine = &SOFT2Matrix::processSampleTable;

    boost::char_separator<char> tdv("\t", "", boost::keep_empty_tokens);
    typedef boost::tokenizer<boost::char_separator<char> > tok_t;
    tok_t tok(aLine, tdv);

    uint32_t n = 0;

    for (tok_t::iterator i = tok.begin(); i != tok.end(); i++, n++)
    {
      if ((*i) == "ID_REF")
        mIdIndex = n;
      else if ((*i) == "VALUE")
        mValueIndex = n;
    }
  }

  void
  processSampleTable(const std::string& aLine)
  {
    if (aLine == "!sample_table_end")
    {
      processLine = &SOFT2Matrix::processSampleIntro;
      return;
    }
  }
};

int
main(int argc, char**argv)
{
  std::string soft, outdir;

  po::options_description desc;

  desc.add_options()
    ("SOFT", po::value<std::string>(&soft), "The SOFT file to process")
    ("outdir", po::value<std::string>(&outdir), "The directory to put the "
     "output into")
    ("help", "produce help message")
    ;

  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("SOFT"))
      wrong = "SOFT";
    else if (!vm.count("outdir"))
      wrong = "outdir";
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;
  if (vm.count("help") || wrong != "")
  {
    std::cout << desc << std::endl;
    return 1;
  }

  if (!fs::is_regular(soft))
  {
    std::cerr << "Invalid SOFT filename supplied." << std::endl;
    return 1;
  }

  if (!fs::is_directory(outdir))
  {
    std::cerr << "Output 'directory' is not a directory." << std::endl;
    return 1;
  }

  SOFT2Matrix s2m(soft, outdir);
  s2m.process();

  return 0;
}
