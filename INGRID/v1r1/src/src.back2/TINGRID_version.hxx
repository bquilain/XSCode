#ifndef TINGRID_version_hxx_seen
#define TINGRID_version_hxx_seen

// Source for TINGRID_version.hxx auto-generated using the
// Tpackage_version.hxx.in template file.

#include <TPackageVersion.hxx>

namespace ND {
    class TINGRID_version;
};

/// This has fields for the library version and can add stuff to the TROOT
/// environment.
class ND::TINGRID_version: public ND::TPackageVersion {
private:
    static TINGRID_version* fThis;
    
public:
    TINGRID_version();
    ~TINGRID_version();

    /// Return a reference to the singleton.
    static TINGRID_version* Get(void); 

    /// Return the version of this library.
    const char* GetName(void) const;

    /// Return the version of this library.
    const char* GetVersion(void) const;

    /// Return the date that this library was compiled.
    const char* GetCompilationDate(void) const;
    
    /// Return the host that this library was compiled on.
    const char* GetCompilationHost(void) const;

    /// Return the directory from which this library was compiled.
    const char* GetCompilationDirectory(void) const;

    /// Return the machine information for the machine that compiled this 
    /// library.  On most machines this is generated from "uname -a".
    const char* GetCompilationMachineInfo(void) const; 

    /// Do any initialization needed for the oaUtility library.  This is
    /// called by the TINGRID_version constructor.
    void Initialize(void);
    
    ClassDef(TINGRID_version,0);
};
#endif
