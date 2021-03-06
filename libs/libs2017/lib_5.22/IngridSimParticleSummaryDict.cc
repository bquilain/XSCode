//
// File generated by rootcint at Fri Apr 14 19:16:02 2017

// Do NOT change. Changes will be lost next time file is generated
//

#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "IngridSimParticleSummaryDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void IngridSimParticleSummary_ShowMembers(void *obj, TMemberInspector &R__insp, char *R__parent);
   static void *new_IngridSimParticleSummary(void *p = 0);
   static void *newArray_IngridSimParticleSummary(Long_t size, void *p);
   static void delete_IngridSimParticleSummary(void *p);
   static void deleteArray_IngridSimParticleSummary(void *p);
   static void destruct_IngridSimParticleSummary(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::IngridSimParticleSummary*)
   {
      ::IngridSimParticleSummary *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::IngridSimParticleSummary >(0);
      static ::ROOT::TGenericClassInfo 
         instance("IngridSimParticleSummary", ::IngridSimParticleSummary::Class_Version(), "../src/IngridSimParticleSummary.h", 9,
                  typeid(::IngridSimParticleSummary), DefineBehavior(ptr, ptr),
                  &::IngridSimParticleSummary::Dictionary, isa_proxy, 0,
                  sizeof(::IngridSimParticleSummary) );
      instance.SetNew(&new_IngridSimParticleSummary);
      instance.SetNewArray(&newArray_IngridSimParticleSummary);
      instance.SetDelete(&delete_IngridSimParticleSummary);
      instance.SetDeleteArray(&deleteArray_IngridSimParticleSummary);
      instance.SetDestructor(&destruct_IngridSimParticleSummary);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::IngridSimParticleSummary*)
   {
      return GenerateInitInstanceLocal((::IngridSimParticleSummary*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::IngridSimParticleSummary*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *IngridSimParticleSummary::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *IngridSimParticleSummary::Class_Name()
{
   return "IngridSimParticleSummary";
}

//______________________________________________________________________________
const char *IngridSimParticleSummary::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IngridSimParticleSummary*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int IngridSimParticleSummary::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IngridSimParticleSummary*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void IngridSimParticleSummary::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IngridSimParticleSummary*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *IngridSimParticleSummary::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IngridSimParticleSummary*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void IngridSimParticleSummary::Streamer(TBuffer &R__b)
{
   // Stream an object of class IngridSimParticleSummary.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> trackid;
      R__b >> parentid;
      R__b >> pdg;
      R__b.ReadStaticArray((float*)momentum);
      R__b.ReadStaticArray((float*)ipos);
      R__b.ReadStaticArray((float*)fpos);
      R__b >> iposflag;
      R__b >> fposflag;
      R__b.ReadStaticArray((float*)dir);
      R__b >> length;
      R__b.CheckByteCount(R__s, R__c, IngridSimParticleSummary::IsA());
   } else {
      R__c = R__b.WriteVersion(IngridSimParticleSummary::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << trackid;
      R__b << parentid;
      R__b << pdg;
      R__b.WriteArray(momentum, 4);
      R__b.WriteArray(ipos, 4);
      R__b.WriteArray(fpos, 4);
      R__b << iposflag;
      R__b << fposflag;
      R__b.WriteArray(dir, 3);
      R__b << length;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void IngridSimParticleSummary::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
      // Inspect the data members of an object of class IngridSimParticleSummary.
      TClass *R__cl = ::IngridSimParticleSummary::IsA();
      Int_t R__ncp = strlen(R__parent);
      if (R__ncp || R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__parent, "trackid", &trackid);
      R__insp.Inspect(R__cl, R__parent, "parentid", &parentid);
      R__insp.Inspect(R__cl, R__parent, "pdg", &pdg);
      R__insp.Inspect(R__cl, R__parent, "momentum[4]", momentum);
      R__insp.Inspect(R__cl, R__parent, "ipos[4]", ipos);
      R__insp.Inspect(R__cl, R__parent, "fpos[4]", fpos);
      R__insp.Inspect(R__cl, R__parent, "iposflag", &iposflag);
      R__insp.Inspect(R__cl, R__parent, "fposflag", &fposflag);
      R__insp.Inspect(R__cl, R__parent, "dir[3]", dir);
      R__insp.Inspect(R__cl, R__parent, "length", &length);
      TObject::ShowMembers(R__insp, R__parent);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_IngridSimParticleSummary(void *p) {
      return  p ? new(p) ::IngridSimParticleSummary : new ::IngridSimParticleSummary;
   }
   static void *newArray_IngridSimParticleSummary(Long_t nElements, void *p) {
      return p ? new(p) ::IngridSimParticleSummary[nElements] : new ::IngridSimParticleSummary[nElements];
   }
   // Wrapper around operator delete
   static void delete_IngridSimParticleSummary(void *p) {
      delete ((::IngridSimParticleSummary*)p);
   }
   static void deleteArray_IngridSimParticleSummary(void *p) {
      delete [] ((::IngridSimParticleSummary*)p);
   }
   static void destruct_IngridSimParticleSummary(void *p) {
      typedef ::IngridSimParticleSummary current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::IngridSimParticleSummary

/********************************************************
* IngridSimParticleSummaryDict.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && (__GNUC__ > 3) && (__GNUC_MINOR__ > 1)
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableIngridSimParticleSummaryDict();

extern "C" void G__set_cpp_environmentIngridSimParticleSummaryDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("../src/IngridSimParticleSummary.h");
  G__cpp_reset_tagtableIngridSimParticleSummaryDict();
}
#include <new>
extern "C" int G__cpp_dllrevIngridSimParticleSummaryDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* IngridSimParticleSummary */
static int G__IngridSimParticleSummaryDict_146_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   IngridSimParticleSummary* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new IngridSimParticleSummary[n];
     } else {
       p = new((void*) gvp) IngridSimParticleSummary[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new IngridSimParticleSummary;
     } else {
       p = new((void*) gvp) IngridSimParticleSummary;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((IngridSimParticleSummary*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) IngridSimParticleSummary::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) IngridSimParticleSummary::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) IngridSimParticleSummary::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      IngridSimParticleSummary::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((IngridSimParticleSummary*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) IngridSimParticleSummary::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) IngridSimParticleSummary::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) IngridSimParticleSummary::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridSimParticleSummaryDict_146_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) IngridSimParticleSummary::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__IngridSimParticleSummaryDict_146_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   IngridSimParticleSummary* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new IngridSimParticleSummary(*(IngridSimParticleSummary*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef IngridSimParticleSummary G__TIngridSimParticleSummary;
static int G__IngridSimParticleSummaryDict_146_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (IngridSimParticleSummary*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((IngridSimParticleSummary*) (soff+(sizeof(IngridSimParticleSummary)*i)))->~G__TIngridSimParticleSummary();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (IngridSimParticleSummary*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((IngridSimParticleSummary*) (soff))->~G__TIngridSimParticleSummary();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__IngridSimParticleSummaryDict_146_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   IngridSimParticleSummary* dest = (IngridSimParticleSummary*) G__getstructoffset();
   *dest = *(IngridSimParticleSummary*) libp->para[0].ref;
   const IngridSimParticleSummary& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* IngridSimParticleSummary */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncIngridSimParticleSummaryDict {
 public:
  G__Sizep2memfuncIngridSimParticleSummaryDict(): p(&G__Sizep2memfuncIngridSimParticleSummaryDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncIngridSimParticleSummaryDict::*p)();
};

size_t G__get_sizep2memfuncIngridSimParticleSummaryDict()
{
  G__Sizep2memfuncIngridSimParticleSummaryDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceIngridSimParticleSummaryDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary))) {
     IngridSimParticleSummary *G__Lderived;
     G__Lderived=(IngridSimParticleSummary*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary),G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableIngridSimParticleSummaryDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<TSchemaHelper>",117,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* IngridSimParticleSummary */
static void G__setup_memvarIngridSimParticleSummary(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary));
   { IngridSimParticleSummary *p; p=(IngridSimParticleSummary*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->trackid)-(long)(p)),105,0,0,-1,-1,-1,1,"trackid=",0,"GEANT4 track ID number");
   G__memvar_setup((void*)((long)(&p->parentid)-(long)(p)),105,0,0,-1,-1,-1,1,"parentid=",0,"GEANT4 parent track ID number");
   G__memvar_setup((void*)((long)(&p->pdg)-(long)(p)),105,0,0,-1,-1,-1,1,"pdg=",0,"PDG particle numbering scheme");
   G__memvar_setup((void*)((long)(&p->momentum)-(long)(p)),102,0,0,-1,-1,-1,1,"momentum[4]=",0,"particle's 4-momentum (GeV)");
   G__memvar_setup((void*)((long)(&p->ipos)-(long)(p)),102,0,0,-1,-1,-1,1,"ipos[4]=",0,"particle's initial position/time (cm, ns)");
   G__memvar_setup((void*)((long)(&p->fpos)-(long)(p)),102,0,0,-1,-1,-1,1,"fpos[4]=",0,"particle's final position/time (cm, ns)");
   G__memvar_setup((void*)((long)(&p->iposflag)-(long)(p)),105,0,0,-1,-1,-1,1,"iposflag=",0,"particle's initial position flag");
   G__memvar_setup((void*)((long)(&p->fposflag)-(long)(p)),105,0,0,-1,-1,-1,1,"fposflag=",0,"particle's final position flag");
   G__memvar_setup((void*)((long)(&p->dir)-(long)(p)),102,0,0,-1,-1,-1,1,"dir[3]=",0,"particle's direction");
   G__memvar_setup((void*)((long)(&p->length)-(long)(p)),102,0,0,-1,-1,-1,1,"length=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarIngridSimParticleSummaryDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncIngridSimParticleSummary(void) {
   /* IngridSimParticleSummary */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary));
   G__memfunc_setup("IngridSimParticleSummary",2472,G__IngridSimParticleSummaryDict_146_0_1, 105, G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Clear",487,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "C - 'Option_t' 10 '\"\"' option", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Print",525,G__IngridSimParticleSummaryDict_146_0_3, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__IngridSimParticleSummaryDict_146_0_4, 85, G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (TClass* (*)())(&IngridSimParticleSummary::Class), 0);
   G__memfunc_setup("Class_Name",982,G__IngridSimParticleSummaryDict_146_0_5, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&IngridSimParticleSummary::Class_Name), 0);
   G__memfunc_setup("Class_Version",1339,G__IngridSimParticleSummaryDict_146_0_6, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (Version_t (*)())(&IngridSimParticleSummary::Class_Version), 0);
   G__memfunc_setup("Dictionary",1046,G__IngridSimParticleSummaryDict_146_0_7, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (void (*)())(&IngridSimParticleSummary::Dictionary), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 2, 1, 1, 0, 
"u 'TMemberInspector' - 1 - insp C - - 0 - parent", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__IngridSimParticleSummaryDict_146_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__IngridSimParticleSummaryDict_146_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&IngridSimParticleSummary::DeclFileName), 0);
   G__memfunc_setup("ImplFileLine",1178,G__IngridSimParticleSummaryDict_146_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (int (*)())(&IngridSimParticleSummary::ImplFileLine), 0);
   G__memfunc_setup("ImplFileName",1171,G__IngridSimParticleSummaryDict_146_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&IngridSimParticleSummary::ImplFileName), 0);
   G__memfunc_setup("DeclFileLine",1152,G__IngridSimParticleSummaryDict_146_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (int (*)())(&IngridSimParticleSummary::DeclFileLine), 0);
   // automatic copy constructor
   G__memfunc_setup("IngridSimParticleSummary", 2472, G__IngridSimParticleSummaryDict_146_0_16, (int) ('i'), 
G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary), -1, 0, 1, 1, 1, 0, "u 'IngridSimParticleSummary' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~IngridSimParticleSummary", 2598, G__IngridSimParticleSummaryDict_146_0_17, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__IngridSimParticleSummaryDict_146_0_18, (int) ('u'), G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary), -1, 1, 1, 1, 1, 0, "u 'IngridSimParticleSummary' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncIngridSimParticleSummaryDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalIngridSimParticleSummaryDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcIngridSimParticleSummaryDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__IngridSimParticleSummaryDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__IngridSimParticleSummaryDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__IngridSimParticleSummaryDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__IngridSimParticleSummaryDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__IngridSimParticleSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__IngridSimParticleSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary = { "IngridSimParticleSummary" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableIngridSimParticleSummaryDict() {
  G__IngridSimParticleSummaryDictLN_TClass.tagnum = -1 ;
  G__IngridSimParticleSummaryDictLN_TBuffer.tagnum = -1 ;
  G__IngridSimParticleSummaryDictLN_TMemberInspector.tagnum = -1 ;
  G__IngridSimParticleSummaryDictLN_TObject.tagnum = -1 ;
  G__IngridSimParticleSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__IngridSimParticleSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableIngridSimParticleSummaryDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__IngridSimParticleSummaryDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__IngridSimParticleSummaryDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__IngridSimParticleSummaryDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__IngridSimParticleSummaryDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__IngridSimParticleSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__IngridSimParticleSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum(&G__IngridSimParticleSummaryDictLN_IngridSimParticleSummary),sizeof(IngridSimParticleSummary),-1,29952,"Simulation (detector mc) particle Summary",G__setup_memvarIngridSimParticleSummary,G__setup_memfuncIngridSimParticleSummary);
}
extern "C" void G__cpp_setupIngridSimParticleSummaryDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupIngridSimParticleSummaryDict()");
  G__set_cpp_environmentIngridSimParticleSummaryDict();
  G__cpp_setup_tagtableIngridSimParticleSummaryDict();

  G__cpp_setup_inheritanceIngridSimParticleSummaryDict();

  G__cpp_setup_typetableIngridSimParticleSummaryDict();

  G__cpp_setup_memvarIngridSimParticleSummaryDict();

  G__cpp_setup_memfuncIngridSimParticleSummaryDict();
  G__cpp_setup_globalIngridSimParticleSummaryDict();
  G__cpp_setup_funcIngridSimParticleSummaryDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncIngridSimParticleSummaryDict();
  return;
}
class G__cpp_setup_initIngridSimParticleSummaryDict {
  public:
    G__cpp_setup_initIngridSimParticleSummaryDict() { G__add_setup_func("IngridSimParticleSummaryDict",(G__incsetup)(&G__cpp_setupIngridSimParticleSummaryDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initIngridSimParticleSummaryDict() { G__remove_setup_func("IngridSimParticleSummaryDict"); }
};
G__cpp_setup_initIngridSimParticleSummaryDict G__cpp_setup_initializerIngridSimParticleSummaryDict;

