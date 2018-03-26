/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program for testing SAMRAI smart pointers
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include <boost/shared_ptr.hpp>

using namespace std;
using namespace SAMRAI;

class A
{
public:
   A() {
   }
   virtual ~A() {
   }
   virtual void
   foo(
      std::string& trace) = 0;
   virtual void
   foo1(
      std::string& trace) = 0;
};

class B
{
public:
   B() {
   }
   virtual ~B() {
   }
   virtual void
   foo(
      std::string& trace) = 0;
   virtual void foo1(
      std::string& trace) {
      trace = "B::foo1()";
   }
};

class C:public A,
   public B
{
public:
   C() {
   }
   virtual ~C() {
   }
   void foo(
      std::string& trace) {
      trace = "C::foo()";
   }
   void foo1(
      std::string& trace) {
      trace = "C::foo1()";
   }
};

class HaveA
{
public:
   HaveA(
      int& fail_count,
      A* a_ptr,
      const A* a_constptr,
      boost::shared_ptr<A> a_smartptr,
      const boost::shared_ptr<A> a_constsmartptr)
   {
      if (a_ptr != a_constptr) {
         fail_count++;
         tbox::perr << "FAILED: - Test #4a: a_ptr != a_constptr" << endl;
      }
      if (a_ptr != a_smartptr.get()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #4b: a_ptr != (A*)a_smartptr" << endl;
      }
      if (a_ptr != a_constsmartptr.get()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #4c: a_ptr != (A*)a_constsmartptr"
                    << endl;
      }
      d_a_ptr = a_ptr;
      d_a_smartptr = a_smartptr;
   }
   virtual ~HaveA() {
   }
   int callFoo()
   {
      int fail_count = 0;
      std::string trace1 = "test-d_a_ptr";
      d_a_ptr->foo(trace1);
      if (trace1 != "C::foo()") {
         fail_count++;
         tbox::perr << "FAILED: - Test #6a: calling d_a_ptr->foo()" << endl;
         tbox::perr << "   Expected C::foo(), Received " << trace1 << endl;
      }
      std::string trace2 = "test-d_a_smartptr";
      d_a_smartptr->foo(trace2);
      if (trace2 != "C::foo()") {
         fail_count++;
         tbox::perr << "FAILED: - Test #6b: calling d_a_smartptr->foo()"
                    << endl;
         tbox::perr << "   Expected C::foo(), Received " << trace2 << endl;
      }
      return fail_count;
   }
   int callFoo1()
   {
      int fail_count = 0;
      string trace = "test-d_a_smartptr";
      d_a_smartptr->foo1(trace);
      if (trace != "C::foo1()") {
         fail_count++;
         tbox::perr << "FAILED: - Test #8a: calling d_a_smartptr->foo1()"
                    << endl;
         tbox::perr << "   Expected C::foo1(), Received " << trace << endl;
      }
      return fail_count;
   }
private:
   A* d_a_ptr;
   boost::shared_ptr<A> d_a_smartptr;
};

class HaveB
{
public:
   HaveB(
      int& fail_count,
      B* b_ptr,
      const B* b_constptr,
      boost::shared_ptr<B> b_smartptr,
      const boost::shared_ptr<B> b_constsmartptr)
   {
      if (b_ptr != b_constptr) {
         fail_count++;
         tbox::perr << "FAILED: - Test #5b: b_ptr != b_constptr" << endl;
      }
      if (b_ptr != b_smartptr.get()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #5b: b_ptr != (A*)b_smartptr" << endl;
      }
      if (b_ptr != b_constsmartptr.get()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #5c: b_ptr != (A*)b_constsmartptr"
                    << endl;
      }
      d_b_ptr = b_ptr;
      d_b_smartptr = b_smartptr;
   }
   virtual ~HaveB() {
   }
   int callFoo()
   {
      int fail_count = 0;
      string trace1 = "test-d_b_ptr";
      d_b_ptr->foo(trace1);
      if (trace1 != "C::foo()") {
         fail_count++;
         tbox::perr << "FAILED: - Test #7a: calling d_b_ptr->foo()" << endl;
         tbox::perr << "   Expected C::foo(), Received " << trace1 << endl;
      }
      std::string trace2 = "test-d_b_smartptr";
      d_b_smartptr->foo(trace2);
      if (trace2 != "C::foo()") {
         fail_count++;
         tbox::perr << "FAILED: - Test #7b: calling d_b_smartptr->foo()"
                    << endl;
         tbox::perr << "   Expected C::foo(), Received " << trace2 << endl;
      }
      return fail_count;
   }
   int callFoo1()
   {
      int fail_count = 0;
      string trace = "test-d_b_smartptr";
      d_b_smartptr->foo1(trace);
      if (trace != "C::foo1()") {
         fail_count++;
         tbox::perr << "FAILED: - Test #9a: calling d_b_smartptr->foo1()"
                    << endl;
         tbox::perr << "   Expected C::foo1(), Received " << trace << endl;
      }
      return fail_count;
   }
private:
   B* d_b_ptr;
   boost::shared_ptr<B> d_b_smartptr;
};

class Derived
{
public:
   Derived() {
   }
   virtual ~Derived() {
   }
};

class ReallyDerived:public Derived
{
public:
   ReallyDerived() {
   }
   virtual ~ReallyDerived() {
   }
};

class ReallyReallyDerived:public ReallyDerived
{
public:
   ReallyReallyDerived() {
   }
   virtual ~ReallyReallyDerived() {
   }
};

template class boost::shared_ptr<A>;
template class boost::shared_ptr<B>;
template class boost::shared_ptr<C>;
template class boost::shared_ptr<Derived>;
template class boost::shared_ptr<ReallyDerived>;
template class boost::shared_ptr<ReallyReallyDerived>;

/*
 * Function to test casting in function call.
 */

int test(
   boost::shared_ptr<ReallyDerived> a,
   boost::shared_ptr<ReallyDerived> b,
   boost::shared_ptr<ReallyDerived> c)
{
   int fail_count = 0;
   if (a) {
      fail_count++;
      tbox::perr << "FAILED: - Test #2a: in test(), a is non-null" << endl;
   }
   if (!b) {
      fail_count++;
      tbox::perr << "FAILED: - Test #2b: in test(), b is null" << endl;
   }
   if (!c) {
      fail_count++;
      tbox::perr << "FAILED: - Test #2c: in test(), c is null" << endl;
   }

   boost::shared_ptr<ReallyReallyDerived> d(
      b,
      boost::detail::dynamic_cast_tag());
   if (d) {
      fail_count++;
      tbox::perr << "FAILED: - Test #2d: in test(), d is non-null" << endl;
   }
   boost::shared_ptr<ReallyReallyDerived> e(
      c,
      boost::detail::dynamic_cast_tag());
   if (!e) {
      fail_count++;
      tbox::perr << "FAILED: - Test #2e: in test(), e is null" << endl;
   }
   return fail_count;
}

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      /*
       * Regular pointer tests.
       */

      boost::shared_ptr<Derived> derived(new Derived());
      boost::shared_ptr<ReallyDerived> really_derived(new ReallyDerived());
      boost::shared_ptr<ReallyReallyDerived> really_really_derived(
         new ReallyReallyDerived());

      /*
       * Test casting to base class Derived.
       */
      boost::shared_ptr<Derived> a(derived);
      if (!a) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1a: a is null" << endl;
      }
      boost::shared_ptr<Derived> b(really_derived);
      if (!b) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1b: b is null" << endl;
      }
      boost::shared_ptr<Derived> c(really_really_derived);
      if (!c) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1c: c is null" << endl;
      }

      /*
       * Test casting to intermediate class ReallyDerived.
       */
      boost::shared_ptr<ReallyDerived> d(
         derived,
         boost::detail::dynamic_cast_tag());
      if (d) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1d: d is non-null" << endl;
      }
      boost::shared_ptr<ReallyDerived> e(really_derived);
      if (!e) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1e: e is null" << endl;
      }
      boost::shared_ptr<ReallyDerived> f(really_really_derived);
      if (!f) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1f: f is null" << endl;
      }

      /*
       * Test casting to most derived class ReallyReallyDerived.
       */
      boost::shared_ptr<ReallyReallyDerived> g(
         derived,
         boost::detail::dynamic_cast_tag());
      if (g) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1g: g is non-null" << endl;
      }
      boost::shared_ptr<ReallyReallyDerived> h(
         really_derived,
         boost::detail::dynamic_cast_tag());
      if (h) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1h: h is non-null" << endl;
      }
      boost::shared_ptr<ReallyReallyDerived> i(really_really_derived);
      if (!i) {
         fail_count++;
         tbox::perr << "FAILED: - Test #1i: i is null" << endl;
      }

      /*
       * Test casting in function call (#2 tests).
       */
      fail_count += test(
         boost::dynamic_pointer_cast<ReallyDerived, Derived>(derived),
         really_derived,
         really_really_derived);

      /*
       * Const pointer tests.
       */

      /*
       * Test casting to base class Derived.
       */
      const boost::shared_ptr<Derived> j(derived);
      if (!j) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3j: j is null" << endl;
      }
      const boost::shared_ptr<Derived> k(really_derived);
      if (!k) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3k: k is null" << endl;
      }
      const boost::shared_ptr<Derived> l(really_really_derived);
      if (!l) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3l: l is null" << endl;
      }

      /*
       * Test casting to intermediate class ReallyDerived.
       */
      const boost::shared_ptr<ReallyDerived> m(
         derived,
         boost::detail::dynamic_cast_tag());
      if (m) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3m: m is non-null" << endl;
      }
      const boost::shared_ptr<ReallyDerived> n(really_derived);
      if (!n) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3n: n is null" << endl;
      }
      const boost::shared_ptr<ReallyDerived> o(really_really_derived);
      if (!o) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3o: o is null" << endl;
      }

      /*
       * Test casting to most derived class ReallyDerived.
       */
      const boost::shared_ptr<ReallyReallyDerived> p(
         derived,
         boost::detail::dynamic_cast_tag());
      if (p) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3p: p is non-null" << endl;
      }
      const boost::shared_ptr<ReallyReallyDerived> q(
         really_derived,
         boost::detail::dynamic_cast_tag());
      if (q) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3q: q is non-null" << endl;
      }
      const boost::shared_ptr<ReallyReallyDerived> r(really_really_derived);
      if (!r) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3r: r is null" << endl;
      }

      /*
       * Test casting const pointer to regular pointers.
       */
#if 0
      const boost::shared_ptr<ReallyReallyDerived> s(derived);
      if (s) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3s: s is non-null" << endl;
      }
      const boost::shared_ptr<ReallyReallyDerived> t(really_derived);
      if (t) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3t: t is non-null" << endl;
      }
      const boost::shared_ptr<ReallyReallyDerived> u(really_really_derived);
      if (u) {
         fail_count++;
         tbox::perr << "FAILED: - Test #3u: u is non-null" << endl;
      }
#endif

      /*
       * Virtual base class pointer tests.
       */

      boost::shared_ptr<C> my_C(new C());

      // #4 tests
      HaveA my_HaveA(fail_count, my_C.get(), my_C.get(), my_C, my_C);

      // #5 tests
      HaveB my_HaveB(fail_count, my_C.get(), my_C.get(), my_C, my_C);

      // #6 tests
      fail_count += my_HaveA.callFoo();

      // #7 tests
      fail_count += my_HaveB.callFoo();

      // #8 tests
      fail_count += my_HaveA.callFoo1();

      // #9 tests
      fail_count += my_HaveB.callFoo1();

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  testpointer" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}
