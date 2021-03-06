//
// File generated by rootcint at Fri Apr 14 19:16:12 2017

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
#include "IngridTrackSummaryDict.h"

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
   void IngridTrackSummary_ShowMembers(void *obj, TMemberInspector &R__insp, char *R__parent);
   static void *new_IngridTrackSummary(void *p = 0);
   static void *newArray_IngridTrackSummary(Long_t size, void *p);
   static void delete_IngridTrackSummary(void *p);
   static void deleteArray_IngridTrackSummary(void *p);
   static void destruct_IngridTrackSummary(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::IngridTrackSummary*)
   {
      ::IngridTrackSummary *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::IngridTrackSummary >(0);
      static ::ROOT::TGenericClassInfo 
         instance("IngridTrackSummary", ::IngridTrackSummary::Class_Version(), "../src/IngridTrackSummary.h", 20,
                  typeid(::IngridTrackSummary), DefineBehavior(ptr, ptr),
                  &::IngridTrackSummary::Dictionary, isa_proxy, 0,
                  sizeof(::IngridTrackSummary) );
      instance.SetNew(&new_IngridTrackSummary);
      instance.SetNewArray(&newArray_IngridTrackSummary);
      instance.SetDelete(&delete_IngridTrackSummary);
      instance.SetDeleteArray(&deleteArray_IngridTrackSummary);
      instance.SetDestructor(&destruct_IngridTrackSummary);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::IngridTrackSummary*)
   {
      return GenerateInitInstanceLocal((::IngridTrackSummary*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::IngridTrackSummary*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *IngridTrackSummary::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *IngridTrackSummary::Class_Name()
{
   return "IngridTrackSummary";
}

//______________________________________________________________________________
const char *IngridTrackSummary::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IngridTrackSummary*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int IngridTrackSummary::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::IngridTrackSummary*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void IngridTrackSummary::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IngridTrackSummary*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *IngridTrackSummary::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::IngridTrackSummary*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void IngridTrackSummary::Streamer(TBuffer &R__b)
{
   // Stream an object of class IngridTrackSummary.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b.ReadStaticArray((float*)vtxi);
      R__b.ReadStaticArray((float*)vtxf);
      R__b >> length;
      R__b >> ekin;
      R__b >> tx;
      R__b >> ty;
      R__b >> etx;
      R__b >> ety;
      R__b >> ex0;
      R__b >> ey0;
      R__b >> covx;
      R__b >> covy;
      R__b >> chi2x;
      R__b >> chi2y;
      R__b >> btheta;
      R__b >> bphi;
      R__b.ReadStaticArray((int*)mrdhitid);
      R__b >> mucl;
      R__b >> vpe;
      R__b >> view;
      R__b >> nhits;
      R__b >> nsimparticles;
      R__b >> nsimemshowers;
      int R__i;
      for (R__i = 0; R__i < 1000; R__i++)
         fIngridHit[R__i].Streamer(R__b);
      for (R__i = 0; R__i < 10; R__i++)
         fSimParticle[R__i].Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, IngridTrackSummary::IsA());
   } else {
      R__c = R__b.WriteVersion(IngridTrackSummary::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b.WriteArray(vtxi, 3);
      R__b.WriteArray(vtxf, 3);
      R__b << length;
      R__b << ekin;
      R__b << tx;
      R__b << ty;
      R__b << etx;
      R__b << ety;
      R__b << ex0;
      R__b << ey0;
      R__b << covx;
      R__b << covy;
      R__b << chi2x;
      R__b << chi2y;
      R__b << btheta;
      R__b << bphi;
      R__b.WriteArray(mrdhitid, 2);
      R__b << mucl;
      R__b << vpe;
      R__b << view;
      R__b << nhits;
      R__b << nsimparticles;
      R__b << nsimemshowers;
      int R__i;
      for (R__i = 0; R__i < 1000; R__i++)
         fIngridHit[R__i].Streamer(R__b);
      for (R__i = 0; R__i < 10; R__i++)
         fSimParticle[R__i].Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void IngridTrackSummary::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
      // Inspect the data members of an object of class IngridTrackSummary.
      TClass *R__cl = ::IngridTrackSummary::IsA();
      Int_t R__ncp = strlen(R__parent);
      if (R__ncp || R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__parent, "vtxi[3]", vtxi);
      R__insp.Inspect(R__cl, R__parent, "vtxf[3]", vtxf);
      R__insp.Inspect(R__cl, R__parent, "length", &length);
      R__insp.Inspect(R__cl, R__parent, "ekin", &ekin);
      R__insp.Inspect(R__cl, R__parent, "tx", &tx);
      R__insp.Inspect(R__cl, R__parent, "ty", &ty);
      R__insp.Inspect(R__cl, R__parent, "etx", &etx);
      R__insp.Inspect(R__cl, R__parent, "ety", &ety);
      R__insp.Inspect(R__cl, R__parent, "ex0", &ex0);
      R__insp.Inspect(R__cl, R__parent, "ey0", &ey0);
      R__insp.Inspect(R__cl, R__parent, "covx", &covx);
      R__insp.Inspect(R__cl, R__parent, "covy", &covy);
      R__insp.Inspect(R__cl, R__parent, "chi2x", &chi2x);
      R__insp.Inspect(R__cl, R__parent, "chi2y", &chi2y);
      R__insp.Inspect(R__cl, R__parent, "btheta", &btheta);
      R__insp.Inspect(R__cl, R__parent, "bphi", &bphi);
      R__insp.Inspect(R__cl, R__parent, "mrdhitid[2]", mrdhitid);
      R__insp.Inspect(R__cl, R__parent, "mucl", &mucl);
      R__insp.Inspect(R__cl, R__parent, "vpe", &vpe);
      R__insp.Inspect(R__cl, R__parent, "view", &view);
      R__insp.Inspect(R__cl, R__parent, "nhits", &nhits);
      R__insp.Inspect(R__cl, R__parent, "nsimparticles", &nsimparticles);
      R__insp.Inspect(R__cl, R__parent, "nsimemshowers", &nsimemshowers);
      R__insp.Inspect(R__cl, R__parent, "fIngridHit[1000]", fIngridHit);
      R__insp.Inspect(R__cl, R__parent, "fSimParticle[10]", fSimParticle);
      TObject::ShowMembers(R__insp, R__parent);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_IngridTrackSummary(void *p) {
      return  p ? new(p) ::IngridTrackSummary : new ::IngridTrackSummary;
   }
   static void *newArray_IngridTrackSummary(Long_t nElements, void *p) {
      return p ? new(p) ::IngridTrackSummary[nElements] : new ::IngridTrackSummary[nElements];
   }
   // Wrapper around operator delete
   static void delete_IngridTrackSummary(void *p) {
      delete ((::IngridTrackSummary*)p);
   }
   static void deleteArray_IngridTrackSummary(void *p) {
      delete [] ((::IngridTrackSummary*)p);
   }
   static void destruct_IngridTrackSummary(void *p) {
      typedef ::IngridTrackSummary current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::IngridTrackSummary

/********************************************************
* IngridTrackSummaryDict.cc
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

extern "C" void G__cpp_reset_tagtableIngridTrackSummaryDict();

extern "C" void G__set_cpp_environmentIngridTrackSummaryDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("../src/IngridTrackSummary.h");
  G__cpp_reset_tagtableIngridTrackSummaryDict();
}
#include <new>
extern "C" int G__cpp_dllrevIngridTrackSummaryDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* IngridTrackSummary */
static int G__IngridTrackSummaryDict_176_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   IngridTrackSummary* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new IngridTrackSummary[n];
     } else {
       p = new((void*) gvp) IngridTrackSummary[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new IngridTrackSummary;
     } else {
       p = new((void*) gvp) IngridTrackSummary;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   IngridTrackSummary* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new IngridTrackSummary(*(IngridTrackSummary*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) IngridTrackSummary(*(IngridTrackSummary*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((IngridTrackSummary*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const IngridTrackSummary*) G__getstructoffset())->NHits());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const IngridTrackSummary*) G__getstructoffset())->NSimParticles());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const IngridTrackSummary*) G__getstructoffset())->NSimEMShowers());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((IngridTrackSummary*) G__getstructoffset())->AddIngridHit((IngridHitSummary*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((const IngridTrackSummary*) G__getstructoffset())->GetIngridHit((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((IngridTrackSummary*) G__getstructoffset())->AddSimParticle((IngridSimParticleSummary*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((const IngridTrackSummary*) G__getstructoffset())->GetSimParticle((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) IngridTrackSummary::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) IngridTrackSummary::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) IngridTrackSummary::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      IngridTrackSummary::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((IngridTrackSummary*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) IngridTrackSummary::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) IngridTrackSummary::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) IngridTrackSummary::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__IngridTrackSummaryDict_176_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) IngridTrackSummary::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef IngridTrackSummary G__TIngridTrackSummary;
static int G__IngridTrackSummaryDict_176_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (IngridTrackSummary*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((IngridTrackSummary*) (soff+(sizeof(IngridTrackSummary)*i)))->~G__TIngridTrackSummary();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (IngridTrackSummary*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((IngridTrackSummary*) (soff))->~G__TIngridTrackSummary();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__IngridTrackSummaryDict_176_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   IngridTrackSummary* dest = (IngridTrackSummary*) G__getstructoffset();
   *dest = *(IngridTrackSummary*) libp->para[0].ref;
   const IngridTrackSummary& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* IngridTrackSummary */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncIngridTrackSummaryDict {
 public:
  G__Sizep2memfuncIngridTrackSummaryDict(): p(&G__Sizep2memfuncIngridTrackSummaryDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncIngridTrackSummaryDict::*p)();
};

size_t G__get_sizep2memfuncIngridTrackSummaryDict()
{
  G__Sizep2memfuncIngridTrackSummaryDict a;
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
extern "C" void G__cpp_setup_inheritanceIngridTrackSummaryDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary))) {
     IngridTrackSummary *G__Lderived;
     G__Lderived=(IngridTrackSummary*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary),G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableIngridTrackSummaryDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<TSchemaHelper>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* IngridTrackSummary */
static void G__setup_memvarIngridTrackSummary(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary));
   { IngridTrackSummary *p; p=(IngridTrackSummary*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->vtxi)-(long)(p)),102,0,0,-1,-1,-1,1,"vtxi[3]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->vtxf)-(long)(p)),102,0,0,-1,-1,-1,1,"vtxf[3]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->length)-(long)(p)),102,0,0,-1,-1,-1,1,"length=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ekin)-(long)(p)),102,0,0,-1,-1,-1,1,"ekin=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->tx)-(long)(p)),102,0,0,-1,-1,-1,1,"tx=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ty)-(long)(p)),102,0,0,-1,-1,-1,1,"ty=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->etx)-(long)(p)),102,0,0,-1,-1,-1,1,"etx=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ety)-(long)(p)),102,0,0,-1,-1,-1,1,"ety=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ex0)-(long)(p)),102,0,0,-1,-1,-1,1,"ex0=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ey0)-(long)(p)),102,0,0,-1,-1,-1,1,"ey0=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->covx)-(long)(p)),102,0,0,-1,-1,-1,1,"covx=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->covy)-(long)(p)),102,0,0,-1,-1,-1,1,"covy=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->chi2x)-(long)(p)),102,0,0,-1,-1,-1,1,"chi2x=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->chi2y)-(long)(p)),102,0,0,-1,-1,-1,1,"chi2y=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->btheta)-(long)(p)),102,0,0,-1,-1,-1,1,"btheta=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->bphi)-(long)(p)),102,0,0,-1,-1,-1,1,"bphi=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->mrdhitid)-(long)(p)),105,0,0,-1,-1,-1,1,"mrdhitid[2]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->mucl)-(long)(p)),102,0,0,-1,-1,-1,1,"mucl=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->vpe)-(long)(p)),102,0,0,-1,-1,-1,1,"vpe=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->view)-(long)(p)),105,0,0,-1,-1,-1,1,"view=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"nhits=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"nsimparticles=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"nsimemshowers=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_TRef),-1,-1,4,"fIngridHit[1000]=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_TRef),-1,-1,4,"fSimParticle[10]=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarIngridTrackSummaryDict() {
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
static void G__setup_memfuncIngridTrackSummary(void) {
   /* IngridTrackSummary */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary));
   G__memfunc_setup("IngridTrackSummary",1856,G__IngridTrackSummaryDict_176_0_1, 105, G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("IngridTrackSummary",1856,G__IngridTrackSummaryDict_176_0_2, 105, G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary), -1, 0, 1, 1, 1, 0, "u 'IngridTrackSummary' - 11 - trk", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Clear",487,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "C - 'Option_t' 10 '\"\"' option", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Print",525,G__IngridTrackSummaryDict_176_0_4, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("NHits",486,G__IngridTrackSummaryDict_176_0_5, 105, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("NSimParticles",1310,G__IngridTrackSummaryDict_176_0_6, 105, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("NSimEMShowers",1268,G__IngridTrackSummaryDict_176_0_7, 105, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AddIngridHit",1163,G__IngridTrackSummaryDict_176_0_8, 121, -1, -1, 0, 1, 1, 1, 0, "U 'IngridHitSummary' - 0 - hit", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetIngridHit",1186,G__IngridTrackSummaryDict_176_0_9, 85, G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridHitSummary), -1, 0, 1, 1, 1, 8, "i - - 0 - i", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("AddSimParticle",1382,G__IngridTrackSummaryDict_176_0_10, 121, -1, -1, 0, 1, 1, 1, 0, "U 'IngridSimParticleSummary' - 0 - part", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetSimParticle",1405,G__IngridTrackSummaryDict_176_0_11, 85, G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridSimParticleSummary), -1, 0, 1, 1, 1, 8, "i - - 0 - i", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__IngridTrackSummaryDict_176_0_12, 85, G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (TClass* (*)())(&IngridTrackSummary::Class), 0);
   G__memfunc_setup("Class_Name",982,G__IngridTrackSummaryDict_176_0_13, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&IngridTrackSummary::Class_Name), 0);
   G__memfunc_setup("Class_Version",1339,G__IngridTrackSummaryDict_176_0_14, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (Version_t (*)())(&IngridTrackSummary::Class_Version), 0);
   G__memfunc_setup("Dictionary",1046,G__IngridTrackSummaryDict_176_0_15, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (void (*)())(&IngridTrackSummary::Dictionary), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 2, 1, 1, 0, 
"u 'TMemberInspector' - 1 - insp C - - 0 - parent", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__IngridTrackSummaryDict_176_0_19, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__IngridTrackSummaryDict_176_0_20, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&IngridTrackSummary::DeclFileName), 0);
   G__memfunc_setup("ImplFileLine",1178,G__IngridTrackSummaryDict_176_0_21, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (int (*)())(&IngridTrackSummary::ImplFileLine), 0);
   G__memfunc_setup("ImplFileName",1171,G__IngridTrackSummaryDict_176_0_22, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&IngridTrackSummary::ImplFileName), 0);
   G__memfunc_setup("DeclFileLine",1152,G__IngridTrackSummaryDict_176_0_23, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (int (*)())(&IngridTrackSummary::DeclFileLine), 0);
   // automatic destructor
   G__memfunc_setup("~IngridTrackSummary", 1982, G__IngridTrackSummaryDict_176_0_24, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__IngridTrackSummaryDict_176_0_25, (int) ('u'), G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary), -1, 1, 1, 1, 1, 0, "u 'IngridTrackSummary' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncIngridTrackSummaryDict() {
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
extern "C" void G__cpp_setup_globalIngridTrackSummaryDict() {
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
}

static void G__cpp_setup_func11() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcIngridTrackSummaryDict() {
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
  G__cpp_setup_func11();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__IngridTrackSummaryDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_TRef = { "TRef" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_IngridHitSummary = { "IngridHitSummary" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_IngridSimParticleSummary = { "IngridSimParticleSummary" , 99 , -1 };
G__linked_taginfo G__IngridTrackSummaryDictLN_IngridTrackSummary = { "IngridTrackSummary" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableIngridTrackSummaryDict() {
  G__IngridTrackSummaryDictLN_TClass.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_TBuffer.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_TMemberInspector.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_TObject.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_TRef.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_IngridHitSummary.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_IngridSimParticleSummary.tagnum = -1 ;
  G__IngridTrackSummaryDictLN_IngridTrackSummary.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableIngridTrackSummaryDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_TRef);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_IngridHitSummary);
   G__get_linked_tagnum_fwd(&G__IngridTrackSummaryDictLN_IngridSimParticleSummary);
   G__tagtable_setup(G__get_linked_tagnum(&G__IngridTrackSummaryDictLN_IngridTrackSummary),sizeof(IngridTrackSummary),-1,30464,"SciBar Track Summary",G__setup_memvarIngridTrackSummary,G__setup_memfuncIngridTrackSummary);
}
extern "C" void G__cpp_setupIngridTrackSummaryDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupIngridTrackSummaryDict()");
  G__set_cpp_environmentIngridTrackSummaryDict();
  G__cpp_setup_tagtableIngridTrackSummaryDict();

  G__cpp_setup_inheritanceIngridTrackSummaryDict();

  G__cpp_setup_typetableIngridTrackSummaryDict();

  G__cpp_setup_memvarIngridTrackSummaryDict();

  G__cpp_setup_memfuncIngridTrackSummaryDict();
  G__cpp_setup_globalIngridTrackSummaryDict();
  G__cpp_setup_funcIngridTrackSummaryDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncIngridTrackSummaryDict();
  return;
}
class G__cpp_setup_initIngridTrackSummaryDict {
  public:
    G__cpp_setup_initIngridTrackSummaryDict() { G__add_setup_func("IngridTrackSummaryDict",(G__incsetup)(&G__cpp_setupIngridTrackSummaryDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initIngridTrackSummaryDict() { G__remove_setup_func("IngridTrackSummaryDict"); }
};
G__cpp_setup_initIngridTrackSummaryDict G__cpp_setup_initializerIngridTrackSummaryDict;

