#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class SOFT2Matrix
{
public:
  SOFT2Matrix(const std::string& aSOFTFile, const std::string& aOutdir)
    : mOutdir(aOutdir), mSOFTFile(aSOFTFile.c_str()), mnSamples(0)
  {
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
  std::vector<std::string> mFields;
  void (SOFT2Matrix::* processLine)(const std::string& aLine);
  uint32_t mnSamples;

  void
  processPlatformIntro(const std::string& aLine)
  {
    if (aLine == "!platform_table_begin")
    {
      processLine = &SOFT2Matrix::processPlatformHeader;
      return;
    }

    if (aLine.substr(0, 22) == "!Platform_sample_id = ")
      mnSamples++;
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
      mFields.push_back(*i);
      if ((*i) == "ID")
        mIdIndex = n;
      else if ((*i) == "Gene Symbol")
        mGeneSymbolIndex = n;

      std::cout << "Platform field: " << (*i) << std::endl;
    }
  }

  uint32_t mIdIndex, mGeneSymbolIndex;

  void
  processPlatformTable(const std::string& aLine)
  {
    if (aLine == "!platform_table_end")
    {
      processLine = &SOFT2Matrix::processSampleIntro;
      return;
    }

    
  }

  void
  processSampleIntro(const std::string& aLine)
  {
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
      std::cout << "Sample field: " << (*i) << std::endl;
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
