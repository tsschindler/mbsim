#include "config.h"
#include <stdlib.h>
#include <iostream>
#include <mbxmlutilshelper/getinstallpath.h>
#include <mbxmlutilshelper/last_write_time.h>
#include <mbxmlutilshelper/dom.h>
#include <xercesc/dom/DOMDocument.hpp>

#include "mbsim/dynamic_system_solver.h"
#include "mbsim/objectfactory.h"
#include "mbsim/integrators/integrator.h"
#include "mbsimflatxml.h"
#include <boost/timer/timer.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string.hpp>

#include <mbsim/analysis/eigenanalysis.h>

#ifndef _WIN32
#  include <dlfcn.h>
#else
#  include <windows.h>
#endif

using namespace std;
using namespace MBXMLUtils;
using namespace xercesc;
using namespace boost;
namespace po = boost::program_options;

namespace {

class SharedLibrary {
  public:
    SharedLibrary(const boost::filesystem::path &file_);
    SharedLibrary(const SharedLibrary& src);
    ~SharedLibrary();
    const boost::filesystem::path file;
    const boost::posix_time::ptime writeTime;
    bool operator<(const SharedLibrary& b) const { return file<b.file; }
  private:
    void init();
#ifndef _WIN32
    void* handle;
#else
    HMODULE handle;
#endif
};

SharedLibrary::SharedLibrary(const boost::filesystem::path &file_) : file(file_),
  writeTime(boost::myfilesystem::last_write_time(file.generic_string())) {
  init();
}

SharedLibrary::SharedLibrary(const SharedLibrary& src) : file(src.file), writeTime(src.writeTime) {
  init();
}

void SharedLibrary::init() {
#ifndef _WIN32
  handle=dlopen(file.generic_string().c_str(), RTLD_NOW | RTLD_LOCAL | RTLD_DEEPBIND);
#else
  handle=LoadLibraryEx(file.generic_string().c_str(), NULL, LOAD_WITH_ALTERED_SEARCH_PATH);
#endif
  if(!handle)
    throw runtime_error("Unable to load the MBSim module: Library '"+file.generic_string()+"' not found.");
}

SharedLibrary::~SharedLibrary() {
#ifndef _WIN32
  dlclose(handle);
#else
  FreeLibrary(handle);
#endif
}

// return the full relative path of a shared library (relative to the install directory, hance including the lib or bin subdir).
// the library base filename 'base' is given with the prefix (lib on Linux) and without the extension.
boost::filesystem::path relLibName(const string &base) {
#ifndef _WIN32
  static const boost::filesystem::path subDir("lib");
  return subDir/("lib"+base+".so.0");
#else
  static const boost::filesystem::path subDir("bin");
  return subDir/("lib"+base+"-0.dll");
#endif
}

// load all MBSim module plugins:
// If a module plugin (shared library) is already loaded but the file has a newer last write time than the
// last write time of the file at the time the shared library was loaded it is unloaded and reloaded.
void loadPlugins() {
  static const NamespaceURI MBSIMPLUGIN("http://mbsim.berlios.de/MBSimPlugin");
  static const boost::filesystem::path installDir(getInstallPath());
  // note: we not not validate the plugin xml files in mbsimflatxml since we do no validated at all in mbsimflatxml (but in mbsimxml)
  static boost::shared_ptr<DOMParser> parser=DOMParser::create(false);

  set<boost::filesystem::path> pluginLibFile;

  // read plugin libraries
  for(boost::filesystem::directory_iterator it=boost::filesystem::directory_iterator(installDir/"share"/"mbsimxml"/"plugins");
      it!=boost::filesystem::directory_iterator(); it++) {
    if(it->path().string().substr(it->path().string().length()-string(".plugin.xml").length())!=".plugin.xml") continue;
    boost::shared_ptr<xercesc::DOMDocument> doc=parser->parse(*it);
    for(xercesc::DOMElement *e=E(E(doc->getDocumentElement())->getFirstElementChildNamed(MBSIMPLUGIN%"libraries"))->
        getFirstElementChildNamed(MBSIMPLUGIN%"Library");
        e!=NULL; e=e->getNextElementSibling())
      pluginLibFile.insert(installDir/relLibName(E(e)->getAttribute("basename")));
  }

  static set<SharedLibrary> loadedPlugin;

  // unload no longer existing plugins or plugins with newer write time
  for(set<SharedLibrary>::iterator it=loadedPlugin.begin(); it!=loadedPlugin.end(); it++)
    if(pluginLibFile.count(it->file)==0 || boost::myfilesystem::last_write_time(it->file.generic_string())>it->writeTime) {
      set<SharedLibrary>::iterator it2=it; it2--;
      loadedPlugin.erase(it);
      it=it2;
    }

  // load plugins which are not already loaded
  for(set<boost::filesystem::path>::iterator it=pluginLibFile.begin(); it!=pluginLibFile.end(); it++)
    loadedPlugin.insert(SharedLibrary(*it));
}

}

namespace MBSim {

int PrefixedStringBuf::sync() {
  // split current buffer into lines
  vector<string> line;
  string full=str();
  boost::split(line, full, boost::is_from_range('\n', '\n'));
  // clear the current buffer
  str("");
  // print each line prefixed with prefix to outstr
  for(vector<string>::iterator it=line.begin(); it!=line.end(); ++it)
    if(it!=--line.end())
      outstr<<prefix<<*it<<"\n";
    else
      if(*it!="")
        outstr<<prefix<<*it;
  outstr<<flush;
  return 0;
}

int MBSimXML::preInit(const po::variables_map &vm, DynamicSystemSolver*& dss, Integrator*& integrator) {

  // setup message streams
  static PrefixedStringBuf infoBuf("Info:    ", cout);
  static PrefixedStringBuf warnBuf("Warning: ", cerr);
  fmatvec::Atom::setCurrentMessageStream(fmatvec::Atom::Info, boost::make_shared<ostream>(&infoBuf));
  fmatvec::Atom::setCurrentMessageStream(fmatvec::Atom::Warn, boost::make_shared<ostream>(&warnBuf));

  loadPlugins();

  // load MBSim project XML document
  shared_ptr<DOMParser> parser=DOMParser::create(false);
  cout << vm["projectfile"].as<string>() << endl;
  shared_ptr<xercesc::DOMDocument> doc=parser->parse(vm["projectfile"].as<string>());
  DOMElement *e=doc->getDocumentElement();

  // create object for DynamicSystemSolver and check correct type
  dss=ObjectFactory::createAndInit<DynamicSystemSolver>(e->getFirstElementChild());

  // create object for Integrator and check correct type
  integrator=ObjectFactory::createAndInit<Integrator>(e->getFirstElementChild()->getNextElementSibling());

  return 0;
}

void MBSimXML::initDynamicSystemSolver(const po::variables_map &vm, DynamicSystemSolver*& dss) {
  if(vm.count("donotintegrate"))
    dss->setTruncateSimulationFiles(false);

  dss->initialize();
}

void MBSimXML::plotInitialState(Integrator* integrator, DynamicSystemSolver* dss) {
  int zSize=dss->getzSize();
  fmatvec::Vec z(zSize);
  if(integrator->getInitialState().size())
    z = integrator->getInitialState();
  else
    dss->initz(z);          
  dss->computeInitialCondition();
  dss->plot(z, 0);
}

void MBSimXML::eigenAnalysis(const po::variables_map &vm, Integrator*& integrator, DynamicSystemSolver*& dss) {
  Eigenanalysis analysis;
//  analysis.setOutputFileName("MBS.eigenanalysis.mat");
  int zSize=dss->getzSize();
  fmatvec::Vec z0(zSize);
  if(integrator->getInitialState().size())
    z0 = integrator->getInitialState();
  else
    dss->initz(z0);          
  analysis.setInitialState(z0);
//  fmatvec::Vec deltaz(dss->getzSize(),fmatvec::INIT,1);
//  analysis.setInitialDeviation(deltaz);
  analysis.setDetermineEquilibriumState(vm.count("equilibrium"));
  analysis.setAutoUpdateAnalysis(vm.count("autoupdate"));
  if(vm.count("amplitude"))
    analysis.setAmplitude(vm["amplitude"].as<double>());
  analysis.setPlotStepSize(integrator->getPlotStepSize());
  string task = vm["task"].as<string>();
  if(task=="eigenanalyis")
    analysis.analyse(*dss);
  if(task=="eigenmode") {
    int i = vm.count("mode")?vm["mode"].as<int>():1;
    analysis.eigenmode(i,*dss);
  }
  else if(task=="eigenmodes")
    analysis.eigenmodes(*dss);
}

void MBSimXML::main(Integrator *integrator, DynamicSystemSolver *dss) {
  boost::timer::cpu_timer t;
  t.start();

  integrator->integrate(*dss);

  t.stop();
  cout<<"Integration CPU times: "<<t.format()<<endl;
}

void MBSimXML::postMain(const po::variables_map &vm, Integrator *&integrator, DynamicSystemSolver*& dss) {

  if(vm.count("savefinalstatevector"))
    dss->writez("statevector.asc", false);
  delete dss;
  delete integrator;
}

}