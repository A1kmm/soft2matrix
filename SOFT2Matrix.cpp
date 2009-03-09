/*
    SOFT2Matrix: Convert from the SOFT format to a packed binary matrix.
    Copyright (C) 2008-2009  Andrew Miller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <list>
#include <boost/tokenizer.hpp>
#include <cstdio>
#include <math.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace bll = boost::lambda;

class SOFT2Matrix
{
public:
  SOFT2Matrix(std::istream& aSOFTFile, const std::string& aOutdir)
    : mOutdir(aOutdir), mSOFTFile(aSOFTFile), mnSamples(0),
      mProbesetCount(0), mProbesets(NULL), mGenes(NULL),
      mGeneProbesetCounts(NULL), mGotSampleTable(true)
  {
    fs::path arrayList(mOutdir);
    arrayList /= "arrays";
    mArrayList = new std::ofstream(arrayList.string().c_str());

    fs::path geneList(mOutdir);
    geneList /= "genes";
    mGeneList = new std::ofstream(geneList.string().c_str());

    fs::path dataFile(mOutdir);
    dataFile /= "data";
    mDataFile = fopen(dataFile.string().c_str(), "w");
  }

  ~SOFT2Matrix()
  {
    if (mProbesets != NULL)
      delete mProbesets;

    if (mGenes != NULL)
      delete mGenes;

    if (mGeneProbesetCounts != NULL)
      delete mGeneProbesetCounts;

    delete mArrayList;
    delete mGeneList;

    if (mDataFile != NULL)
      fclose(mDataFile);
  }

  void
  process()
  {
    std::string l;
    processLine = &SOFT2Matrix::processPlatformIntro;

    while (mSOFTFile.good())
    {
      std::getline(mSOFTFile, l);
      (this->*processLine)(l);
    }

    if (mNextId != mSampleIds.end())
    {
      std::cout << "There were samples indicated in the platform sample " 
                << "list but missing in the data file."<< std::endl;
    }
  }

  void
  loadHGNCDatabase(const std::string& aPath)
  {
    io::filtering_istream db;
    db.push(io::file_source(aPath));

    // Skip the header...
    std::string entry;
    std::getline(db, entry);

    while (db.good())
    {
      std::getline(db, entry);

      boost::tokenizer<boost::char_separator<char> >
        tok(entry, boost::char_separator<char>("\t", "",
                                               boost::keep_empty_tokens));
      std::vector<std::string> v(tok.begin(), tok.end());

      if (v.size() < 6)
        continue;

      if (v[3] != "Approved")
        continue;

      uint32_t hgncId = strtoul(v[0].c_str(), NULL, 10);
      addHGNCMapping(v[1], hgncId, true);

      static const boost::regex rtok("[, ]+");

      addHGNCMapping(v[2], hgncId, false);

      boost::sregex_token_iterator rti1
        (make_regex_token_iterator(v[4], rtok, -1));
      boost::sregex_token_iterator end;
      for (; rti1 != end; rti1++)
        addHGNCMapping(*rti1, hgncId, false);

      boost::sregex_token_iterator rti2
        (make_regex_token_iterator(v[5], rtok, -1));
      for (; rti2 != end; rti2++)
        addHGNCMapping(*rti2, hgncId, false);
    }
  }

private:
  fs::path mOutdir;
  std::istream& mSOFTFile;
  std::ofstream *mArrayList, *mGeneList;
  FILE * mDataFile;
  void (SOFT2Matrix::* processLine)(const std::string& aLine);
  uint32_t mnSamples;
  std::list<std::string> mSampleIds;
  std::list<std::string>::iterator mNextId;
  double* mProbesets, * mGenes;
  uint32_t* mGeneProbesetCounts;
  bool mGotSampleTable;

  void
  processPlatformIntro(const std::string& aLine)
  {
    if (aLine == "!platform_table_begin")
    {
      mNextId = mSampleIds.begin();
      processLine = &SOFT2Matrix::processPlatformHeader;
      return;
    }

    if (aLine.substr(0, 22) == "!Platform_sample_id = ")
    {
      std::string sampleId(aLine.substr(22));
      mSampleIds.push_back(sampleId);
      (*mArrayList) << sampleId << std::endl;
      mnSamples++;
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
  uint32_t mProbesetCount, mGeneCount;
  std::map<uint32_t, uint32_t> mGeneIndexByHGNCId;

  void
  fillProbesetArrayWithNans()
  {
    for (uint32_t i = 0; i < mProbesetCount; i++)
      mProbesets[i] = std::numeric_limits<double>::quiet_NaN();
  }

  uint32_t
  findHGNCIdByName(const std::string& aName, bool stripDashes = true)
  {
    // Look up the name from HGNC...
    std::map<std::string, uint32_t>::iterator i
      (mHGNCIdMappings.find(aName));
    if (i != mHGNCIdMappings.end())
      return (*i).second;

    // See if it ends in a number...
    static const boost::regex endNumber("(\\-?)([0-9]+)$");
    boost::smatch res;
    if (boost::regex_search(aName, res, endNumber))
    {
      i = mHGNCIdMappings.find(res.prefix().str());
      if (i != mHGNCIdMappings.end())
        return (*i).second;

      std::string tryAlso;
      if (res[2].str() == "alpha")
        tryAlso = "A";
      else if (res[2].str() == "beta")
        tryAlso = "B";
      else if (res[2].str() == "1")
        tryAlso = "I";
      else if (res[2].str() == "2")
        tryAlso = "II";

      std::string attempt(res.prefix().str());
      attempt += tryAlso;
      i = mHGNCIdMappings.find(attempt);
      if (i != mHGNCIdMappings.end())
        return (*i).second;
    }

    // Try adding a suffix like 1 or A...
    std::string attempt = aName + "1";
    i = mHGNCIdMappings.find(attempt);
    if (i != mHGNCIdMappings.end())
      return (*i).second;
    
    attempt = aName + "A";
    i = mHGNCIdMappings.find(attempt);
    if (i != mHGNCIdMappings.end())
      return (*i).second;

    if (stripDashes)
    {
      // Strip out all dashes and repeat...
      std::string dashless(boost::replace_all_copy(aName, "ALPHA", "A"));
      boost::replace_all(dashless, "-", "");
      return findHGNCIdByName(dashless, false);
    }

    return 0;
  }

  void
  platformTableDone()
  {
    mProbesets = new double[mProbesetCount];
    fillProbesetArrayWithNans();
    
    mGeneCount = mUsedHGNCIds.size();
    mGenes = new double[mGeneCount];
    mGeneProbesetCounts = new uint32_t[mGeneCount];

    std::cout << "mGeneCount = " << mGeneCount << std::endl
              << "mnSamples = " << mnSamples << std::endl;

    std::map<uint32_t, uint32_t> hgncIdToGeneIndex;
    uint32_t geneIndex(0);

    typedef std::pair<uint32_t, uint32_t> pairu32;

    for (std::set<uint32_t>::iterator i = mUsedHGNCIds.begin();
         i != mUsedHGNCIds.end();
         i++, geneIndex++)
    {
      hgncIdToGeneIndex.insert(pairu32(*i, geneIndex));
      (*mGeneList) << mNameByHGNCId[*i] << std::endl;
    }

    std::transform(
                   mProbesetHGNCIdList.begin(),
                   mProbesetHGNCIdList.end(),
                   std::back_inserter(mProbesetGeneList),
                   bll::bind<pairu32>
                   (
                    bll::constructor<pairu32>(),
                    bll::bind(&pairu32::first, bll::_1),
                    bll::var(hgncIdToGeneIndex)
                    [bll::bind<uint32_t>(&pairu32::second, bll::_1)]
                   )
                  );

    processLine = &SOFT2Matrix::processSampleIntro;
  }

  void
  processPlatformTable(const std::string& aLine)
  {
    if (aLine == "!platform_table_end")
    {
      platformTableDone();
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

    if (symbol == "")
      return;

    // Symbol is a list of genes, some of which will be in HGNC...
    static const boost::regex geneSep(" // ");
    boost::sregex_token_iterator rti
      (boost::make_regex_token_iterator(symbol, geneSep, -1));

    std::set<uint32_t> seenIds;
    for (; rti != boost::sregex_token_iterator(); rti++)
    {
      uint32_t hgncid(findHGNCIdByName(*rti));
      if (hgncid == 0)
        continue;

      if (seenIds.count(hgncid) != 0)
        continue;
      seenIds.insert(hgncid);
      mUsedHGNCIds.insert(hgncid);

      mProbesetHGNCIdList.push_back(std::pair<uint32_t, uint32_t>
                                    (hgncid, mProbesetCount));
    }

    mProbesetIndexById.insert(std::pair<std::string, uint32_t>(id, mProbesetCount));
    mProbesetCount++;
  }

  std::map<std::string, uint32_t> mProbesetIndexById;
  std::list<std::pair<uint32_t, uint32_t> > mProbesetHGNCIdList, mProbesetGeneList;
  std::set<uint32_t> mUsedHGNCIds;

  void
  processSampleIntro(const std::string& aLine)
  {
    if (aLine.substr(0, 10) == "^SAMPLE = ")
    {
      if (!mGotSampleTable)
      {
        // This means we found two ^SAMPLE records with no intervening 
        // !sample_table_begin lines! Write a message...
        std::cout << "Warning: Next sample found without a "
          "!sample_table_begin line!" << std::endl;

        // Next, we need to write out a NaN-filled placeholder entry for the
        // missing data...
        for (uint32_t i = 0; i < mGeneCount; i++)
          mGenes[i] = std::numeric_limits<double>::quiet_NaN();
        if (fwrite(mGenes, mGeneCount * sizeof(double), 1, mDataFile) != 1)
        {
          std::cout << "Failed to write record." << std::endl;
        }

        fflush(mDataFile); // Makes checking file sizes easier...
      }

      mGotSampleTable = false;
      std::string sampId = aLine.substr(10);
      if (sampId != *mNextId)
        std::cout << "Sample ID mismatch: expected "
                  << *mNextId << " got " << sampId
                  << std::endl;
      else
        std::cout << "Proc: " << sampId << std::endl;
      mNextId++;
      return;
    }
    if (aLine == "!sample_table_begin")
    {
      mGotSampleTable = true;
      processLine = &SOFT2Matrix::processSampleHeader;
    }
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
  sampleTableDone()
  {
    memset(mGenes, 0, sizeof(double) * mGeneCount);
    memset(mGeneProbesetCounts, 0, sizeof(uint32_t) * mGeneCount);

    for (std::list<std::pair<uint32_t, uint32_t> >::iterator i(mProbesetGeneList.begin());
         i != mProbesetGeneList.end();
         i++)
    {
      uint32_t probeset = (*i).first;
      uint32_t gene = (*i).second;

      if (isfinite(mProbesets[probeset]))
      {
        mGeneProbesetCounts[gene]++;
        mGenes[gene] += mProbesets[probeset];
      }
    }
    
    for (uint32_t i = 0; i < mGeneCount; i++)
    {
      if (mGeneProbesetCounts[i] == 0)
        mGenes[i] = std::numeric_limits<double>::quiet_NaN();
      else
        mGenes[i] /= mGeneProbesetCounts[i];
    }

    if (fwrite(mGenes, mGeneCount * sizeof(double), 1, mDataFile) != 1)
    {
      std::cout << "Failed to write record." << std::endl;
    }

    fflush(mDataFile); // Makes checking file sizes easier...

    fillProbesetArrayWithNans();
    processLine = &SOFT2Matrix::processSampleIntro;
  }

  void
  processSampleTable(const std::string& aLine)
  {
    if (aLine == "!sample_table_end")
    {
      sampleTableDone();
      return;
    }

    boost::char_separator<char> tdv("\t", "", boost::keep_empty_tokens);
    typedef boost::tokenizer<boost::char_separator<char> > tok_t;
    tok_t tok(aLine, tdv);

    uint32_t n = 0;

    std::string id, value;

    for (tok_t::iterator i = tok.begin(); i != tok.end(); i++, n++)
    {
      if (n == mIdIndex)
        id = *i;
      else if (n == mValueIndex)
        value = *i;
    }
    
    if (mProbesetIndexById.count(id) == 0)
    {
      // std::cout << "Unknown probe ID " << id << std::endl;
      return;
    }

    mProbesets[mProbesetIndexById[id]] = strtod(value.c_str(), NULL);
  }

  std::string
  cleanup_HGNC_name(const std::string& aName)
  {
    std::string uc(boost::algorithm::to_upper_copy(aName));
    boost::algorithm::replace_all(uc, "-", "");

    return uc;
  }

  std::map<std::string, uint32_t> mHGNCIdMappings;
  std::map<uint32_t, std::string> mNameByHGNCId;

  void addHGNCMapping(const std::string& aMapping, uint32_t aHGNC,
                      bool aOverride)
  {
    std::string dcmapping(cleanup_HGNC_name(aMapping));

    if (aOverride)
      mNameByHGNCId.insert(std::pair<uint32_t, std::string>(aHGNC, aMapping));

    std::map<std::string, uint32_t>::iterator i =
      mHGNCIdMappings.find(dcmapping);
    if (i != mHGNCIdMappings.end())
    {
      if (!aOverride)
        return;

      mHGNCIdMappings.erase(i);
    }

    mHGNCIdMappings.insert(std::pair<std::string, uint32_t>
                           (dcmapping, aHGNC));
  }
};

int
main(int argc, char**argv)
{
  std::string soft, outdir, hgnc;

  po::options_description desc;

  desc.add_options()
    ("SOFT", po::value<std::string>(&soft), "The SOFT file to process")
    ("outdir", po::value<std::string>(&outdir), "The directory to put the "
     "output into")
    ("hgnc", po::value<std::string>(&hgnc), "File containing the HGNC names database")
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
    else if (!vm.count("hgnc"))
      wrong = "hgnc";
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

  if (!fs::is_regular(hgnc))
  {
    std::cerr << "Invalid HGNC filename supplied." << std::endl;
    return 1;
  }

  io::filtering_istream str;
  str.push(io::bzip2_decompressor());
  str.push(io::file_source(soft));

  {
    SOFT2Matrix s2m(str, outdir);
    s2m.loadHGNCDatabase(hgnc);
    s2m.process();
  }

  return 0;
}
