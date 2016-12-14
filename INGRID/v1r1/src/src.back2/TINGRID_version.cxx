#include "TOADatabase.hxx"

// Source for TINGRID_version.cxx auto-generated using the
// Tpackage_version.cxx.in template file.

#include "TINGRID_version.hxx"
#include "INGRID_version.h"

ClassImp(ND::TINGRID_version);

// Trickiness so that the package version is automatically added to the
// list of used packages.
static ND::TINGRID_version INGRID_version;

ND::TINGRID_version* ND::TINGRID_version::fThis = NULL;

ND::TINGRID_version::TINGRID_version() {
    fThis = ND::TINGRID_version::Get();
}

ND::TINGRID_version::~TINGRID_version() {}

void ND::TINGRID_version::Initialize(void) {
    // register this package.
    ND::TOADatabase::Get().PackageSet().insert(fThis);
}

ND::TINGRID_version* ND::TINGRID_version::Get(void) {
    // Make sure that fThis is initialized;
    if (!fThis) {
        // Make sure that fThis is not null before allocating a real pointer.
        // This cruft is required so that there isn't an infinite recursion
        // while fThis is initialized.
        fThis = (ND::TINGRID_version*) 1;
        // Allocate real space for the fThis pointer.
        fThis = new ND::TINGRID_version;
        // Now initialize
        fThis->Initialize();
    }
    // Return the pointer.
    return fThis;
}

const char* ND::TINGRID_version::GetName(void) const {
    return INGRID_NAME;
}

const char* ND::TINGRID_version::GetVersion(void) const {
    return INGRID_VERSION;
}

const char* ND::TINGRID_version::GetCompilationDate(void) const {
    return INGRID_COMPILE_DATE;
}

const char* ND::TINGRID_version::GetCompilationHost(void) const {
    return INGRID_COMPILE_HOST;
}

const char* ND::TINGRID_version::GetCompilationDirectory(void) const {
    return INGRID_COMPILE_DIR;
}

const char* ND::TINGRID_version::GetCompilationMachineInfo(void) const {
    return INGRID_COMPILE_UNAME;
}
